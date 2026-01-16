import json
import logging
import os
import subprocess
import sys
import warnings
from collections import defaultdict
from copy import deepcopy
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

import pandas as pd
from cobra.io import load_json_model, load_matlab_model, read_sbml_model
from platformdirs import user_data_path

from . import data
from .conversion import (
    ConvAgora,
    ConvBase,
    ConvBigg,
    ConvCarveme,
    ConvGapseq,
    ConvMetanetx,
    ConvModelseed,
    no_changes_for_notconv,
    remove_zero_for_notconv,
    replace_square_brackets,
)
from .creation import SuperModel
from .curation import get_duplicated_reactions, remove_b_type_exchange
from .dbs import download_db, get_bigg_network
from .genes import (
    get_final_fasta_with_ncbi_assemble,
    get_genes_gapseq,
    get_genes_not_gapseq,
)
from .periplasmic import getSuggestionPeriplasmic
from .selection import run_selection
from .structural import runStructuralConversion, runSuggestionsMet


def add_charge_mass_info(df_mbs):
    with open(files(data) / "masses_and_charges.json") as fh:
        df_bigg_extra = pd.DataFrame(
            [
                (bigg_id, *x)
                for bigg_id, info in json.load(fh).items()
                if info is not None
                for x in zip(info["formula"], info["mass"], info["charges"])
                if x[1]
                is not None  # some formulas are not valid i.e. contain R, X or Z
            ],
            columns=["universal_bigg_id", "formula", "mass", "charge"],
        ).drop_duplicates(subset="universal_bigg_id", keep="first")
    return df_mbs.merge(df_bigg_extra, on="universal_bigg_id", how="left").assign(
        formula=lambda x: x.formula.fillna("")
    )


def get_env():
    """
    Function to appen additional path to system PATH env var and return dict
    with new env. In case the package is used in conda env.
    """

    sep = ":"
    if sys.platform == "win32":
        sep = ";"

    add_bin = str(Path(sys.executable).parent)
    tmp_env = os.environ.copy()
    tmp_env["PATH"] = f"{tmp_env['PATH']}{sep}{add_bin}"
    return tmp_env


class LoggerContext:
    def __init__(self, logger_name, show_logs=False):
        self.__show_logs__ = show_logs
        self.__logger__ = logging.getLogger(logger_name)

    def __enter__(self):
        if not self.__show_logs__:
            self.__logger__.setLevel(logging.CRITICAL)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        if not self.__show_logs__:
            self.__logger__.setLevel(logging.NOTSET)


@lru_cache
def load_sbml_model(path_to_model, show_logs: bool = False):

    # Read the cobra model
    with LoggerContext("cobra", show_logs):
        if str(path_to_model).endswith(".json"):
            model = load_json_model(path_to_model)
        elif str(path_to_model).endswith(".mat"):
            model = load_matlab_model(path_to_model)
        else:
            model = read_sbml_model(path_to_model)

    return model


class GatheredModels:
    """
    Class that gathers information and necessary conversion results for all
    models. Input for the class and tool in general is dictionary with all
    models and related information.

    Parameters
    ----------
    dict_of_all_models_with_feature : dict
        Dictionary of the following format:
        {
            model_id:
            {
                'path_to_model': str,
                'model_type': str one of (agora, carveme, gapseq, modelseed) or custom,
                'path_to_genome': str (can be '' or None if convert_genes = False)
            }
        }
        Note: if custom `model_type` create type class in advance,
    assembly : None, optional
    path_final_genome_nt : None, optional
    path_final_genome_aa : None, optional
    custom_model_type : None, optional

    Notes
    -----
    If all 3 parameters bellow are None then gene conversion is not done and
    genomes for model_id are not need.
    """

    def __init__(
        self,
        custom_model_type=None,
        clear_db_cache=False,
    ):
        # If specified, clear the cached conversion tables and dictionaries
        if clear_db_cache:
            for p in user_data_path("gemsembler").iterdir():
                p.unlink()

        self.__conf = {
            "agora": {
                "remove_b": False,
                "db_name": "weird_bigg",
                "wo_periplasmic": True,
                "conv_strategy": ConvAgora(),
                "genome_model_strategy": get_genes_not_gapseq,
                "alter_notconv_m": replace_square_brackets,
                "alter_notconv_r": no_changes_for_notconv,
            },
            "carveme": {
                "remove_b": False,
                "db_name": "bigg",
                "wo_periplasmic": False,
                "conv_strategy": ConvCarveme(),
                "genome_model_strategy": get_genes_not_gapseq,
                "alter_notconv_m": no_changes_for_notconv,
                "alter_notconv_r": no_changes_for_notconv,
            },
            "bigg": {
                "remove_b": False,
                "db_name": "bigg",
                "wo_periplasmic": False,
                "conv_strategy": ConvBigg(),
                "genome_model_strategy": get_genes_not_gapseq,
                "alter_notconv_m": no_changes_for_notconv,
                "alter_notconv_r": no_changes_for_notconv,
            },
            "metanetx": {
                "remove_b": False,
                "db_name": "bigg",
                "wo_periplasmic": False,  # TODO: need to check this
                "conv_strategy": ConvMetanetx(),
                "genome_model_strategy": get_genes_not_gapseq,
                "alter_notconv_m": no_changes_for_notconv,
                "alter_notconv_r": no_changes_for_notconv,
            },
            "gapseq": {
                "remove_b": False,
                "db_name": "modelseed",
                "wo_periplasmic": False,
                "conv_strategy": ConvGapseq(),
                "genome_model_strategy": get_genes_gapseq,
                "alter_notconv_m": remove_zero_for_notconv,
                "alter_notconv_r": remove_zero_for_notconv,
            },
            "modelseed": {
                "remove_b": True,
                "db_name": "modelseed",
                "wo_periplasmic": True,
                "conv_strategy": ConvModelseed(),
                "genome_model_strategy": get_genes_not_gapseq,
                "alter_notconv_m": remove_zero_for_notconv,
                "alter_notconv_r": remove_zero_for_notconv,
            },
        }

        self.__models = {}
        self.first_stage_selected_metabolites = None
        self.first_stage_selected_reactions = None
        self.structural_first_run_reactions = defaultdict(dict)
        self.second_stage_selected_reactions = None
        self.structural_first_run_metabolites = defaultdict(dict)
        self.second_stage_selected_metabolites = None
        self.structural_second_run_reactions = defaultdict(dict)
        self.third_stage_selected_reactions = None
        self.many_to_one_sug = defaultdict(dict)
        self.periplasmic_metabolites = defaultdict(dict)
        self.periplasmic_reactions = defaultdict(dict)

    def __contains__(self, item):
        return item in self.__models

    def __len__(self):
        return len(self.__models)

    def get_conf(self, model_type=None):
        if model_type is None:
            return deepcopy(self.__conf)
        else:
            return deepcopy(self.__conf.get(model_type))

    def get_model_attrs(self, model_id=None, attr=None):
        if model_id is not None and attr is not None:
            return deepcopy(self.__models.get(model_id).get(attr))
        elif model_id is not None and attr is None:
            return deepcopy(self.__models.get(model_id))
        else:
            return deepcopy(self.__models)

    @property
    def same_db_models(self):
        same_db_models = defaultdict(dict)
        for model_id, model_attrs in self.__models.items():
            model_type = model_attrs["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            same_db_models[db_name][model_id] = model_type
        return same_db_models

    @property
    def converted_metabolites(self):
        conv_mbs = defaultdict(dict)
        for model_id, model_attrs in self.__models.items():
            converter = self.__conf.get(model_attrs["model_type"]).get("conv_strategy")
            conv_mbs[model_id] = converter.convert_all_metabolites(
                model_attrs["preprocess_model"]
            )
        return conv_mbs

    @property
    def converted_reactions(self):
        conv_rcts = defaultdict(dict)
        for model_id, model_attrs in self.__models.items():
            converter = self.__conf.get(model_attrs["model_type"]).get("conv_strategy")
            conv_rcts[model_id] = converter.convert_all_reactions(
                model_attrs["preprocess_model"]
            )
        return conv_rcts

    def run(self):
        # run first convertion
        print("Running initial convertion")
        self.first_stage_selected_metabolites = run_selection(
            self.same_db_models, self.converted_metabolites, "highest"
        )
        self.first_stage_selected_reactions = run_selection(
            self.same_db_models, self.converted_reactions, "highest"
        )

        # run first structural conversion
        print("Running 1st structural convertion")
        bigg_network = get_bigg_network()
        for model_id, first_sel in self.first_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            self.structural_first_run_reactions[model_id] = runStructuralConversion(
                db_name,
                first_sel,
                self.first_stage_selected_metabolites[model_id],
                self.__models[model_id]["preprocess_model"],
                bigg_network,
                False,
            )
        # run second stage selection for first structural reactions
        self.second_stage_selected_reactions = run_selection(
            self.same_db_models,
            self.structural_first_run_reactions,
            "structural",
            replace_with_consistent=False,
        )

        # get suggestions from structural reactions for metabolites
        print("Running structural suggestions for metabolites")
        for model_id, rs_struct_sel in self.second_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_mod = self.__conf.get(model_type).get("db_name")
            (
                self.structural_first_run_metabolites[model_id],
                self.many_to_one_sug[model_id],
            ) = runSuggestionsMet(
                db_mod,
                self.structural_first_run_reactions[model_id],
                rs_struct_sel,
                self.first_stage_selected_metabolites[model_id],
            )
        # run second stage selection for suggestions for metabolites from structural
        self.second_stage_selected_metabolites = run_selection(
            self.same_db_models, self.structural_first_run_metabolites, "structural"
        )

        # run second structural conversion with suggestions for metabolites
        print("Running 2d structural convertion")
        for model_id, sec_sel in self.second_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            self.structural_second_run_reactions[model_id] = runStructuralConversion(
                db_name,
                sec_sel,
                self.second_stage_selected_metabolites[model_id],
                self.__models[model_id]["preprocess_model"],
                bigg_network,
                self.__conf.get(model_type).get("wo_periplasmic"),
            )
        # run third stage selection for first structural reactions
        self.third_stage_selected_reactions = run_selection(
            self.same_db_models,
            self.structural_second_run_reactions,
            "structural",
            replace_with_consistent=False,
        )

        print("Introducing periplasmic compartment")
        # introducing periplasmic compartment for models, that don't have it originally
        for model_id, th_sel in self.third_stage_selected_reactions.items():
            if self.__conf.get(self.__models[model_id]["model_type"]).get(
                "wo_periplasmic"
            ):
                (
                    self.periplasmic_metabolites[model_id],
                    self.periplasmic_reactions[model_id],
                ) = getSuggestionPeriplasmic(
                    th_sel,
                    self.structural_second_run_reactions[model_id],
                    self.second_stage_selected_metabolites[model_id],
                    self.__models[model_id]["preprocess_model"],
                    bigg_network,
                )
            else:
                (
                    self.periplasmic_metabolites[model_id],
                    self.periplasmic_reactions[model_id],
                ) = ({}, {})

    def get_input_dictionaries(self):
        final_r_sel = defaultdict(dict)
        final_r_not_sel = defaultdict(dict)
        final_m_sel = defaultdict(dict)
        final_m_not_sel = defaultdict(dict)
        periplasmic_m = defaultdict(dict)
        periplasmic_r = defaultdict(dict)
        for model_id in self.__models.keys():
            final_r_sel[model_id] = {}
            final_r_not_sel[model_id] = {}
            final_m_sel[model_id] = {}
            final_m_not_sel[model_id] = {}
            for orig_r_id, sel_r in self.third_stage_selected_reactions[
                model_id
            ].items():
                if (sel_r.to_one_id == True and sel_r.from_one_id == True) or (
                    sel_r.to_one_id == True
                    and sel_r.from_one_id == False
                    and (
                        (
                            orig_r_id
                            in self.__models[model_id]["duplicated_reactions"][
                                "ID"
                            ].tolist()
                        )
                        or (sel_r.highest_consistent == ["Biomass"])
                    )
                ):
                    final_r_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, sel_r.highest_consistent]}
                    )
                else:
                    orig_r_id_alt = self.__conf[self.__models[model_id]["model_type"]][
                        "alter_notconv_r"
                    ](orig_r_id)
                    final_r_not_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, [orig_r_id_alt]]}
                    )
                if (
                    len(
                        self.__models[model_id]["preprocess_model"]
                        .reactions.get_by_id(orig_r_id)
                        .reactants
                    )
                    > 24
                ):
                    final_r_not_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, sel_r.highest_consistent]}
                    )
            for orig_m_id, sel_m in self.second_stage_selected_metabolites[
                model_id
            ].items():
                if orig_m_id in self.periplasmic_metabolites[model_id]:
                    if self.periplasmic_metabolites[model_id][orig_m_id].replace:
                        final_m_sel[model_id].update(
                            {
                                orig_m_id: [
                                    ["p"],
                                    [
                                        self.periplasmic_metabolites[model_id][
                                            orig_m_id
                                        ].bigg_p
                                    ],
                                ]
                            }
                        )
                    else:
                        final_m_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, sel_m.highest_consistent]}
                        )
                        periplasmic_m[model_id].update(
                            {
                                orig_m_id: [
                                    ["p"],
                                    [
                                        self.periplasmic_metabolites[model_id][
                                            orig_m_id
                                        ].bigg_p
                                    ],
                                ]
                            }
                        )
                else:
                    if sel_m.to_one_id == True and sel_m.from_one_id == True:
                        final_m_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, sel_m.highest_consistent]}
                        )
                    else:
                        orig_m_id_alt = self.__conf[
                            self.__models[model_id]["model_type"]
                        ]["alter_notconv_m"](orig_m_id)
                        final_m_not_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, [orig_m_id_alt]]}
                        )
            if self.periplasmic_reactions[model_id]:
                for p_r_id, p_r in self.periplasmic_reactions[model_id].items():
                    periplasmic_r[model_id].update({p_r_id: p_r.metabolites_changed_p})

        return (
            final_m_sel,
            final_m_not_sel,
            final_r_sel,
            final_r_not_sel,
            periplasmic_m,
            periplasmic_r,
        )

    def assemble_supermodel(
        self,
        output_folder: str,
        assembly_id=None,
        path_final_genome_nt=None,
        path_final_genome_aa=None,
        evalue_threshold=0.001,
        do_old_genes=None,
        do_mix_conv_notconv=False,
        and_as_solid=False,
        do_old_locus_tag=True,
    ):
        # Check if assembly and final genome are present.
        # If not, throw a warning.
        if do_old_genes is None:
            do_old_genes = {model_id: False for model_id in self.__models.keys()}
        if not assembly_id and not path_final_genome_nt and not path_final_genome_aa:
            warnings.warn(
                "\nWarning! No final genome for gene conversion is provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide either assembly id or custom "
                "fasta files (nt/aa/both), \nto which genes must be converted."
            )
            gene_path = None
        elif assembly_id and (path_final_genome_nt or path_final_genome_aa):
            warnings.warn(
                "\nWarning! Both assembly and user final genome for gene conversion are provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide one of both either assembly id or custom "
                "fasta files (nt/aa/both), to which genes must be converted."
            )
            gene_path = None
            if do_old_genes is None:
                do_old_genes = {model_id: True for model_id in self.__models.keys()}
        else:
            # GGE: Check that BLAST is installed
            proc = subprocess.run(
                f"makeblastdb -h",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if proc.returncode != 0:
                raise OSError(
                    "Check that makeblastdb (and BLAST in general) is properly installed!"
                    "\n(a remainder, that you NEED to install blast to use this package)"
                )
            output_folder = Path(output_folder)
            gene_path = output_folder / "tmp_gene_conversion"
            gene_path.mkdir(exist_ok=True, parents=True)
            db_path = gene_path / "blast_db"
            db_path.mkdir(exist_ok=True, parents=True)
            print("Downloading assembly from NCBI")
            if assembly_id:
                (
                    path_final_genome_nt,
                    path_final_genome_aa,
                ) = get_final_fasta_with_ncbi_assemble(
                    output_folder, assembly_id, do_old_locus_tag=do_old_locus_tag
                )
            if path_final_genome_nt is not None:
                print("Building BLAST database")
                run = subprocess.run(
                    f"makeblastdb -in {path_final_genome_nt} -out "
                    f"{Path(db_path, 'nt_db')} -dbtype nucl "
                    "-title nt_db -parse_seqids",  # WindowsFix ('' removed)  [seems to work on Linux as well]
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=get_env(),
                )
                if run.returncode != 0:
                    raise OSError(f"Failed to run makeblastdb: {run.stderr.decode()}")
            if path_final_genome_aa is not None:
                print("Building BLAST database")
                run = subprocess.run(
                    f"makeblastdb -in {path_final_genome_aa} -out "
                    f"{Path(db_path, 'aa_db')} -dbtype prot "
                    "-title aa_db -parse_seqids",  # WindowsFix ('' removed) [seems to work on Linux as well]
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=get_env(),
                )
                if run.returncode != 0:
                    raise OSError(f"Failed to run makeblastdb: {run.stderr.decode()}")
            for model_id, model_data in self.__models.items():
                path_to_genome = model_data["path_to_genome"]
                if path_to_genome is None or path_to_genome == "":
                    continue
                print(f"Running gene conversion with BLAST for {model_id}")
                out_blast_file = gene_path / (model_id + "_blast.tsv")
                model_gene_file, aa_status = self.__conf[model_data["model_type"]][
                    "genome_model_strategy"
                ](
                    gene_path,
                    model_data["path_to_genome"],
                    model_data["preprocess_model"],
                    model_data["model_type"],
                    model_id,
                )
                blast_command = ""
                db_name = ""
                if aa_status and path_final_genome_aa is not None:
                    blast_command = "blastp"
                    db_name = "aa_db"  # WindowsFix ('' removed) [seems to work on Linux as well]
                elif aa_status and path_final_genome_nt is not None:
                    blast_command = "tblastn"
                    db_name = "nt_db"  # WindowsFix ('' removed) [seems to work on Linux as well]
                elif not aa_status and path_final_genome_nt is not None:
                    blast_command = "blastn"
                    db_name = "nt_db"  # WindowsFix ('' removed) [seems to work on Linux as well]
                elif not aa_status and path_final_genome_aa is not None:
                    blast_command = "blastx"
                    db_name = "aa_db"  # WindowsFix ('' removed) [seems to work on Linux as well]
                if blast_command == "" or db_name == "":
                    warnings.warn("\nWarning! Something wrong with aa/nt in files/DB")
                elif not model_gene_file:
                    warnings.warn("\nWarning! Something wrong with gene file")
                else:
                    run = subprocess.run(
                        f"{blast_command} -query {model_gene_file} "
                        f"-db {Path(db_path, db_name)} "
                        f"-max_target_seqs 1 -evalue {evalue_threshold} "
                        f"-outfmt 6 -out {out_blast_file}",  # WindowsFix ('' removed) [seems to work on Linux as well]
                        shell=True,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        env=get_env(),
                    )
                    if run.returncode != 0:
                        raise OSError(
                            f"Failed to run {blast_command}: {run.stderr.decode()}"
                        )
        print("Assembling Supermodel")
        # Get final tables to create new objects
        (
            final_m_sel,
            final_m_not_sel,
            final_r_sel,
            final_r_not_sel,
            periplasmic_m,
            periplasmic_r,
        ) = self.get_input_dictionaries()

        # Create supermodel
        bigg_data_m = download_db(
            "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
            "bigg_models_metabolites.txt.gz",
        ).pipe(add_charge_mass_info)

        bigg_data_r = download_db(
            "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt",
            "bigg_models_reactions.txt.gz",
        )

        supermodel = SuperModel(
            False,
            {
                "type": "SuperModel",
                "args": {
                    "final_m_sel": final_m_sel,
                    "final_m_not_sel": final_m_not_sel,
                    "final_r_sel": final_r_sel,
                    "final_r_not_sel": final_r_not_sel,
                    "all_models_data": self.__models,
                    "additional_periplasmic_m": periplasmic_m,
                    "periplasmic_r": periplasmic_r,
                    "m_db_info": bigg_data_m,
                    "r_db_info": bigg_data_r,
                    "gene_folder": gene_path,
                    "do_old_genes": do_old_genes,
                    "do_mix_conv_notconv": do_mix_conv_notconv,
                    "and_as_solid": and_as_solid,
                },
            },
        )
        return supermodel

    def set_configuration(
        self,
        model_type: str,
        remove_b: bool,
        db_name: str,
        wo_periplasmic: bool,
        conv_strategy,
        genome_model_strategy,
        **kwargs,
    ):
        assert isinstance(conv_strategy, ConvBase)

        # TODO: add checks on conf input arg
        self.__conf[model_type] = {
            "remove_b": remove_b,
            "db_name": db_name,
            "wo_periplasmic": wo_periplasmic,
            "conv_strategy": conv_strategy,
            "genome_model_strategy": genome_model_strategy,
            **kwargs,
        }

    def add_model(
        self,
        model_id: str,
        path_to_model: str,
        model_type: str,
        path_to_genome: str = None,
        show_logs: bool = False,
    ):
        # Run checks on model_id and model_type
        # TODO: check with conversion dictionaries
        assert model_id not in self.__models, f"model_id {model_id} already used"
        assert model_type in self.__conf, f"Missing configuration for {model_type}"

        model = load_sbml_model(path_to_model, show_logs)

        # Populate the internal data
        self.__models[model_id] = {
            "original_model": deepcopy(model),
            "path_to_model": path_to_model,
            "model_type": model_type,
            "path_to_genome": path_to_genome,
        }

        # If model_type requires it, remove `_b` extensions
        if self.__conf.get(model_type).get("remove_b"):
            model = remove_b_type_exchange(model)
        self.__models[model_id]["preprocess_model"] = model

        dupl_r = get_duplicated_reactions(model)
        self.__models[model_id]["duplicated_reactions"] = dupl_r

    def add_models_and_run(self, models_list):
        for model in models_list:
            self.add_model(**model)
        self.run()
