import itertools
import operator
import re
import resource
import sys
import warnings
from collections import defaultdict
from math import ceil
from os.path import exists
from pathlib import PosixPath

import dill
import pandas as pd

from .comparison import (
    getCore,
    getCoreCoefficients,
    getCoreConnections,
    getCoreGPR,
    getCoreLowerBounds,
    getCoreUpperBounds,
    getDifference,
)
from .genes import makeNewGPR, uniteGPR


class KnowledgeConnectingOldNew:
    """Gathering methods to connect old and new models"""

    def __init__(
        self,
        m_dictionaries: list,
        r_dictionaries: list,
        m_periplasmic: dict,
        r_periplasmic: dict,
        gene_folder: PosixPath,
    ):
        self.m_go_old_new = m_dictionaries[0]
        self.m_go_new_old = m_dictionaries[1]
        self.m_go_old_new_nc = m_dictionaries[2]
        self.m_go_new_old_nc = m_dictionaries[3]
        self.r_go_old_new = r_dictionaries[0]
        self.r_go_new_old = r_dictionaries[1]
        self.r_go_old_new_nc = r_dictionaries[2]
        self.r_go_new_old_nc = r_dictionaries[3]
        self.m_periplasmic = m_periplasmic
        self.r_periplasmic = r_periplasmic
        self.periplasmic_models = list(r_periplasmic.keys())
        self.g_conversion_tables = defaultdict()
        for model_id in self.m_go_old_new.keys():
            if gene_folder is None:
                self.g_conversion_tables[model_id] = None
                continue
            blast_file = gene_folder / (model_id + "_blast.tsv")
            try:
                conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            except:
                self.g_conversion_tables[model_id] = None
            else:
                conversion_table.columns = [
                    "old_id",
                    "new_id",
                    "identity",
                    "length",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8",
                    "9",
                    "10",
                    "11",
                ]
                self.g_conversion_tables[model_id] = conversion_table

    def get_old_mets(self, model_id: str, new_id: str, do_notconv: bool):
        if do_notconv:
            return self.m_go_new_old_nc.get(new_id).get(model_id)
        else:
            return self.m_go_new_old.get(new_id).get(model_id)

    def get_new_mets(self, model_id: str, old_id: str, do_notconv: bool):
        if do_notconv:
            return self.m_go_old_new_nc.get(model_id).get(old_id)
        else:
            return self.m_go_old_new.get(model_id).get(old_id)

    def get_p_met_ids(self, model_id: str):
        return list(self.m_periplasmic.get(model_id, {}).keys())

    def get_old_rs(self, model_id: str, old_id: str, do_notconv: bool):
        if do_notconv:
            return self.r_go_new_old_nc.get(old_id).get(model_id)
        else:
            return self.r_go_new_old.get(old_id).get(model_id)

    def get_new_rs(self, model_id: str, old_id: str, do_notconv: bool):
        if do_notconv:
            return self.r_go_old_new_nc.get(model_id).get(old_id)
        else:
            return self.r_go_old_new.get(model_id).get(old_id)

    def get_p_rs(self, model_id: str):
        return list(self.r_periplasmic.get(model_id, {}).keys())

    def get_p_mets_for_p_r(self, model_id: str, old_id_r: str):
        return list(self.r_periplasmic.get(model_id).get(old_id_r, {}).keys())

    def get_new_gene_id(self, model_id: str, old_id: str):
        if self.g_conversion_tables[model_id] is None:
            return old_id
        attr_new = self.g_conversion_tables[model_id][
            self.g_conversion_tables[model_id]["old_id"] == old_id
        ]["new_id"]
        if not attr_new.empty:
            return attr_new.values[0]
        else:
            return "not_found"


class NewElement:
    """ New object class - one metabolite or reaction for supermodel. """

    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
    ):
        self.id = new_id
        self.compartments = {"assembly": compartments}
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        self.converted = converted
        for ps in possible_sources:
            if ps == source:
                self.compartments.update({ps: compartments})
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.compartments.update({ps: []})
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def _update_new_element(
        self, id_to_update: str, compart_to_update: [str], source: str,
    ):
        self.sources.update({source: self.sources.get(source) + 1})
        if source not in self.in_models["models_list"]:
            self.in_models["models_amount"] = self.in_models["models_amount"] + 1
            self.in_models["models_list"].append(source)
        self.annotation.get(source).append(id_to_update)
        self.compartments.update(
            {source: self.compartments.get(source) + compart_to_update}
        )
        for c in compart_to_update:
            if c not in self.compartments["assembly"]:
                self.compartments["assembly"].append(c)


class NewMetabolite(NewElement):
    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
        m_database_info: pd.core.frame.DataFrame,
    ):
        super().__init__(
            new_id, old_id, compartments, source, possible_sources, converted
        )
        if converted:
            id_noc = re.sub("_([cep])$", "", new_id)
            name = m_database_info[m_database_info["universal_bigg_id"] == id_noc][
                "name"
            ].values[0]
        else:
            name = "Not converted"
        self.name = name
        self.reactions = {k: [] for k in possible_sources}
        self.reactions.update({"assembly": [], "comparison": {}})
        self.formula = {k: [] for k in possible_sources}

    def _find_reactions(self, connections: KnowledgeConnectingOldNew, do_notconv=False):
        for model_id in self.sources.keys():
            old_mets = connections.get_old_mets(model_id, self.id, do_notconv)
            if not old_mets:
                continue
            new_r = []
            for old_met in old_mets:
                p_met_ids = connections.get_p_met_ids(model_id)
                if old_met.id in p_met_ids:
                    # Old metabolite has additional periplasmic version in supermodel.
                    # So we need to split its reactions between them
                    for reaction in old_met.reactions:
                        new_rs_1old = connections.get_new_rs(
                            model_id, reaction.id, do_notconv
                        )
                        if not new_rs_1old:
                            continue
                        met_is_periplasmic = self.id.endswith("_p")
                        old_p_mets_for_p_r = connections.get_p_mets_for_p_r(
                            model_id, reaction.id
                        )
                        met_was_converted_to_periplasmic = (
                            old_met.id in old_p_mets_for_p_r
                        )
                        if (met_is_periplasmic & met_was_converted_to_periplasmic) or (
                            (not met_is_periplasmic)
                            & (not met_was_converted_to_periplasmic)
                        ):
                            new_r.append(new_rs_1old[0])
                else:
                    # No need to split reactions
                    for r in old_met.reactions:
                        new_rs_1old = connections.get_new_rs(model_id, r.id, do_notconv)
                        if new_rs_1old:
                            new_r.append(new_rs_1old[0])
            if new_r:
                self.reactions[model_id] = list(set(new_r))


class NewReaction(NewElement):
    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
        r_database_info: pd.core.frame.DataFrame,
    ):
        super().__init__(
            new_id, old_id, compartments, source, possible_sources, converted
        )
        if converted:
            id_noc = new_id.replace("sink_", "DM_")
            name = r_database_info[r_database_info["bigg_id"] == id_noc]["name"]
            if (not name.empty) and (not name.isnull().values.any()):
                name = name.values[0]
            else:
                name = ""
            equation = r_database_info[r_database_info["bigg_id"] == id_noc][
                "reaction_string"
            ]
            if not equation.empty:
                equation = equation.values[0]
            else:
                equation = None
        else:
            name = "Not converted"
            equation = None
        self.name = name
        self.reaction = equation
        base_keys = possible_sources + ["assembly"]
        self.reactants = {k: [] for k in base_keys}
        self.reactants.update({"comparison": {}})
        self.products = {k: [] for k in base_keys}
        self.products.update({"comparison": {}})
        self.metabolites = {k: {} for k in base_keys + ["comparison"]}
        self.lower_bound = {k: [] for k in base_keys}
        self.lower_bound.update({"comparison": {}})
        self.upper_bound = {k: [] for k in base_keys}
        self.upper_bound.update({"comparison": {}})
        self.subsystem = {k: [] for k in possible_sources}
        self.genes = {k: [] for k in base_keys}
        self.genes.update({"comparison": {}})
        self.gene_reaction_rule = {k: [] for k in base_keys}
        self.gene_reaction_rule.update({k + "_mixed": [] for k in possible_sources})
        self.gene_reaction_rule.update({"comparison": {}})

    def __sel_met_from_p_model_for_p_r(
        self,
        old_react_metabolites: list,
        old_react_id: str,
        model_id: str,
        connections: KnowledgeConnectingOldNew,
        do_notconv: bool,
    ):
        out_met = {}
        for met in old_react_metabolites:
            new_mets = connections.get_new_mets(model_id, met.id, do_notconv)
            if not new_mets:
                continue
            if len(new_mets) == 1:
                out_met.update({new_mets[0]: met})
            else:
                test_of_single_entry = True
                for new_met in new_mets:
                    new_met_is_periplasmic = new_met.id.endswith("_p")
                    p_mets_for_p_r = connections.get_p_mets_for_p_r(
                        model_id, old_react_id
                    )
                    new_met_was_converted_to_periplasmic = met.id in p_mets_for_p_r
                    if (
                        new_met_is_periplasmic & new_met_was_converted_to_periplasmic
                    ) or (
                        (not new_met_is_periplasmic)
                        & (not new_met_was_converted_to_periplasmic)
                    ):
                        out_met.update({new_met: met})
                        if test_of_single_entry:
                            test_of_single_entry = False
                            continue
                        if not test_of_single_entry:
                            problem_m = [n.id for n, m in out_met.items() if m == met]
                            warnings.warn(
                                f"Something went wrong with periplasmic"
                                f"connections of {self.id}. Problematic"
                                f"metabolites are {' '.join(problem_m)}"
                            )
        return out_met

    def __sel_met_from_p_model_for_not_p_r(
        self,
        old_react_metabolites: list,
        model_id: str,
        connections: KnowledgeConnectingOldNew,
        do_notconv: bool,
    ):
        out_met = {}
        for met in old_react_metabolites:
            new_mets = connections.get_new_mets(model_id, met.id, do_notconv)
            if not new_mets:
                continue
            if len(new_mets) == 1:
                out_met.update({new_mets[0]: met})
            else:
                test_of_single_entry = True
                for new_met in new_mets:
                    if not new_met.id.endswith("_p"):
                        out_met.update({new_met: met})
                        if test_of_single_entry:
                            test_of_single_entry = False
                            continue
                        if not test_of_single_entry:
                            problem_m = [n.id for n, m in out_met.items() if m == met]
                            warnings.warn(
                                f"Something went wrong with periplasmic"
                                f"connections of {self.id}. Problematic"
                                f"metabolites are {' '.join(problem_m)}"
                            )
        return out_met

    def _find_reactants_products(
        self, connections: KnowledgeConnectingOldNew, m_type: str, do_notconv=False,
    ):
        for model_id in self.sources.keys():
            old_react = connections.get_old_rs(model_id, self.id, do_notconv)
            # old_react is list, usually with 1 element,
            # but even if not 1, these reactions have the same r equation,
            # so we can take any and I tool the 1st
            if not old_react:
                continue
            old_react_react_prod = getattr(old_react[0], m_type)
            model_has_periplasmic_changes = model_id in connections.periplasmic_models
            p_rs = connections.get_p_rs(model_id)
            reaction_has_periplasmic_changes = old_react[0].id in p_rs
            if (model_has_periplasmic_changes) & (reaction_has_periplasmic_changes):
                met_for_p_r = self.__sel_met_from_p_model_for_p_r(
                    old_react_react_prod,
                    old_react[0].id,
                    model_id,
                    connections,
                    do_notconv,
                )
                for m in met_for_p_r.keys():
                    getattr(self, m_type).get(model_id).append(m)
            elif model_has_periplasmic_changes & (not reaction_has_periplasmic_changes):
                met_for_not_p_r = self.__sel_met_from_p_model_for_not_p_r(
                    old_react_react_prod, model_id, connections, do_notconv,
                )
                for m in met_for_not_p_r.keys():
                    getattr(self, m_type).get(model_id).append(m)
            else:
                # There was no periplasmic perturbation in the model
                # Only 1 element in new_reacts_prods is expected
                for react_prod in old_react_react_prod:
                    new_reacts_prods = connections.get_new_mets(
                        model_id, react_prod.id, do_notconv
                    )
                    if new_reacts_prods:
                        getattr(self, m_type).get(model_id).append(new_reacts_prods[0])
                        if len(new_reacts_prods) > 1:
                            warnings.warn(
                                f"Unexpected not unique connections between new "
                                f"and old metabolite without periplasmic story."
                                f"Model: {model_id}. Old metabolite: {react_prod.id}."
                                f"New metabolites: "
                                f"{' '.join([n.id for n in new_reacts_prods])}."
                                f"But {new_reacts_prods[0]} was selected."
                            )

    def _find_metabolites(
        self, connections: KnowledgeConnectingOldNew, do_notconv=False
    ):
        for model_id in self.sources.keys():
            old_react = connections.get_old_rs(model_id, self.id, do_notconv)
            # old_react is list, usually with 1 element,
            # but even if not 1, these reactions have the same r equation,
            # so we can take any and I tool the 1st
            if not old_react:
                continue
            old_react_metabolites = old_react[0].metabolites
            model_has_periplasmic_changes = model_id in connections.periplasmic_models
            p_rs = connections.get_p_rs(model_id)
            reaction_has_periplasmic_changes = old_react[0].id in p_rs
            if model_has_periplasmic_changes & reaction_has_periplasmic_changes:
                met_for_p_r = self.__sel_met_from_p_model_for_p_r(
                    list(old_react_metabolites.keys()),
                    old_react[0].id,
                    model_id,
                    connections,
                    do_notconv,
                )
                for m, v in met_for_p_r.items():
                    self.metabolites.get(model_id).update({m: old_react_metabolites[v]})
            elif model_has_periplasmic_changes & (not reaction_has_periplasmic_changes):
                met_for_not_p_r = self.__sel_met_from_p_model_for_not_p_r(
                    list(old_react_metabolites.keys()),
                    model_id,
                    connections,
                    do_notconv,
                )
                for m, v in met_for_not_p_r.items():
                    self.metabolites.get(model_id).update({m: old_react_metabolites[v]})
            else:
                # There was no periplasmic perturbation in the model
                # Only 1 element in new_mets is expected
                for met, koef in old_react_metabolites.items():
                    new_mets = connections.get_new_mets(model_id, met.id, do_notconv)
                    if new_mets:
                        self.metabolites.get(model_id).update({new_mets[0]: koef})
                        if len(new_mets) > 1:
                            warnings.warn(
                                f"Unexpected not unique connections between new "
                                f"and old metabolite without periplasmic story."
                                f"Model: {model_id}. Old metabolite: {met.id}."
                                f"New metabolites: {' '.join([n.id for n in new_mets])}"
                            )

    def _find_gene_and_gpr(
        self, connections: KnowledgeConnectingOldNew, gene_folder, do_notconv=False
    ):
        genes_to_add = defaultdict(list)
        for model_id in self.in_models["models_list"]:
            genes_not_to_add = ["not_found"]
            old_rs = connections.get_old_rs(model_id, self.id, do_notconv)
            if not old_rs:
                continue
            new_gpr_unite_r = []
            for oldr in old_rs:
                if not oldr.genes:
                    continue
                gene_convert = {}
                for oldrg in oldr.genes:
                    pot_new_g_id = connections.get_new_gene_id(model_id, oldrg.id)
                    gene_convert.update({oldrg.id: pot_new_g_id})
                    if (not do_notconv) & (
                        connections.g_conversion_tables[model_id] is not None
                    ):
                        genes_not_to_add.append(oldrg.id)
                    if pot_new_g_id not in genes_not_to_add:
                        genes_to_add[model_id].append(pot_new_g_id)
                        genes_not_to_add.append(pot_new_g_id)
                old_gpr = oldr.gene_reaction_rule
                new_gpr, mix_gpr = makeNewGPR(old_gpr, gene_convert)
                if new_gpr:
                    new_gpr_unite_r.append(new_gpr)
                self.gene_reaction_rule.get(model_id + "_mixed").append(mix_gpr)
            if len(new_gpr_unite_r) == 1:
                self.gene_reaction_rule.get(model_id).append(new_gpr_unite_r[0])
            elif len(new_gpr_unite_r) >= 1:
                united_gpr = uniteGPR(new_gpr_unite_r)
                self.gene_reaction_rule.get(model_id).append(united_gpr)
        return genes_to_add


class SetofNewElements:
    """ Setting dictionaries for all metabolites or reactions:
    selected for supermodel - self.assembly and not selected - self.notconverted. """

    def __add_new_elements(
        self,
        element_type: str,
        selected: dict,
        where_to_add: str,
        model_ids: list,
        convered: bool,
        db_info: pd.core.frame.DataFrame,
    ):
        new_elements = {"metabolites": NewMetabolite, "reactions": NewReaction}
        for mod_id in model_ids:
            if mod_id not in list(selected.keys()):
                continue
            objects = selected.get(mod_id)
            for key in objects.keys():
                for new_id in objects[key][1]:
                    comp = objects[key][0]
                    if new_id in getattr(self, where_to_add).keys():
                        if convered == getattr(self, where_to_add)[new_id].converted:
                            getattr(self, where_to_add).get(new_id)._update_new_element(
                                key, comp, mod_id
                            )
                        elif (
                            new_id + "_convert_" + str(convered)
                            in getattr(self, where_to_add).keys()
                        ):
                            getattr(self, where_to_add).get(
                                new_id + "_convert_" + str(convered)
                            )._update_new_element(key, comp, mod_id)
                        else:
                            new_conv = new_elements[element_type](
                                new_id + "_convert_" + str(convered),
                                key,
                                comp,
                                mod_id,
                                model_ids,
                                convered,
                                db_info,
                            )
                            getattr(self, where_to_add).update(
                                {new_id + "_convert_" + str(convered): new_conv}
                            )
                    else:
                        new = new_elements[element_type](
                            new_id, key, comp, mod_id, model_ids, convered, db_info
                        )
                        getattr(self, where_to_add).update({new_id: new})

    def __init__(
        self,
        element_type: str,
        selected: dict,
        not_selected: dict,
        model_ids: [str],
        db_info: pd.core.frame.DataFrame,
        do_mix_conv_notconv: bool,
        additional=None,
    ):
        self.assembly = {}
        for source in selected.keys():
            setattr(self, source, {})
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        self.__add_new_elements(
            element_type, selected, "assembly", model_ids, True, db_info
        )
        if additional:
            self.__add_new_elements(
                element_type, additional, "assembly", model_ids, True, db_info
            )
        if do_mix_conv_notconv:
            self.__add_new_elements(
                element_type, not_selected, "assembly", model_ids, False, db_info
            )
        else:
            self.__add_new_elements(
                element_type, not_selected, "notconverted", model_ids, False, db_info
            )
        for new_id, new_obj in self.assembly.items():
            for model_id in new_obj.in_models["models_list"]:
                getattr(self, model_id).update({new_id: new_obj})

    def _makeForwardBackward(
        self,
        all_models: dict,
        selected: dict,
        obj_type: "metabolites" or "reactions",
        additional=None,
        not_selected=None,
    ):
        """ Creating dictionaries linking metabolites/reactions:
            NewObject in supermodel with old original ID and OldObject in original models with new ID in supermodel """
        go_old_new = defaultdict(dict)
        go_new_old = defaultdict(dict)
        go_old_new_notconv = defaultdict(dict)
        go_new_old_notconv = defaultdict(dict)
        model_ids = list(selected.keys())
        for model_id in model_ids:
            go_old_new[model_id] = {}
            go_old_new_notconv[model_id] = {}
            for key, value in selected.get(model_id).items():
                new_obj = [self.assembly.get(value[1][0])]
                if additional:
                    if key in list(additional.get(model_id, {}).keys()):
                        new_obj.append(
                            self.assembly.get(additional.get(model_id).get(key)[1][0])
                        )
                go_old_new[model_id].update({key: new_obj})
            if not_selected is not None:
                for key, value in not_selected.get(model_id).items():
                    new_obj_nc = [self.notconverted.get(value[1][0])]
                    go_old_new_notconv[model_id].update({key: new_obj_nc})
        sel_not_sel = {"assembly": go_new_old}
        if not_selected is not None:
            sel_not_sel.update({"notconverted": go_new_old_notconv})
        for s, d_to_add in sel_not_sel.items():
            for k, v in getattr(self, s).items():
                for mod_id in v.in_models["models_list"]:
                    old_ids = v.annotation[mod_id]
                    old_obj = [
                        getattr(
                            all_models[mod_id]["preprocess_model"], obj_type
                        ).get_by_id(i)
                        for i in old_ids
                    ]
                    d_to_add[k].update({mod_id: old_obj})
        return [go_old_new, go_new_old, go_old_new_notconv, go_new_old_notconv]


class NewGene(object):
    """Class for one gene with new or old locus tag as ID and IDs from original models in annotation"""

    def __init__(
        self,
        new_id: str,
        old_id: str,
        source: str,
        possible_sources: [str],
        converted: bool,
    ):
        self.id = new_id
        self.sources = {}
        self.converted = converted
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        self.reactions = {"assembly": [], "comparison": {}}
        for ps in possible_sources:
            self.reactions.update({ps: []})
            if ps == source:
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def _updateNewGene(self, id_to_update: str, source: str):
        self.sources.update({source: self.sources.get(source) + 1})
        if source not in self.in_models["models_list"]:
            self.in_models["models_amount"] = self.in_models["models_amount"] + 1
            self.in_models["models_list"].append(source)
        self.annotation.get(source).append(id_to_update)

    def _find_reactions(
        self,
        all_models_data: dict,
        connections: KnowledgeConnectingOldNew,
        do_conv=False,
    ):
        for model_id in self.in_models["models_list"]:
            old_g_ids = self.annotation.get(model_id)
            for old_g_id in old_g_ids:
                oldg_r_ids = [
                    gr.id
                    for gr in all_models_data[model_id]["preprocess_model"]
                    .genes.get_by_id(old_g_id)
                    .reactions
                ]
                for r_id in oldg_r_ids:
                    new_rs = connections.get_new_rs(model_id, r_id, do_conv)
                    if new_rs:
                        for new_r in new_rs:
                            if new_r not in self.reactions.get(model_id):
                                self.reactions.get(model_id).append(new_r)


class SetofNewGenes(object):
    """ Setting dictionaries for all genes selected for supermodel - self.converted and not selected - self.notconverted. """

    def __addNewGenes_conv(
        self, all_models_data: dict, gene_folder: PosixPath, do_mix_conv_notconv: bool
    ):
        for model_id in list(all_models_data.keys()):
            blast_file = gene_folder / (model_id + "_blast.tsv")
            try:
                conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            except:
                warnings.warn(
                    f"\nWarning! File {str(blast_file)} can't be opened."
                    f"\nOld gene will be used"
                )
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                        getattr(self, model_id).update(
                            {gene.id: self.assembly.get(gene.id)}
                        )
                    else:
                        new_gene = NewGene(
                            gene.id,
                            gene.id,
                            model_id,
                            list(all_models_data.keys()),
                            False,
                        )
                        self.assembly.update({gene.id: new_gene})
                        getattr(self, model_id).update({gene.id: new_gene})
            else:
                conversion_table.columns = [
                    "old_id",
                    "new_id",
                    "identity",
                    "length",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8",
                    "9",
                    "10",
                    "11",
                ]
                if do_mix_conv_notconv:
                    to_add = "assembly"
                else:
                    to_add = "notconverted"
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    old_gene_id = gene.id
                    attr = conversion_table[conversion_table["old_id"] == old_gene_id][
                        "new_id"
                    ]
                    if attr.empty:
                        if gene.id in getattr(self, to_add).keys():
                            getattr(self, to_add).get(gene.id)._updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                gene.id,
                                gene.id,
                                model_id,
                                list(all_models_data.keys()),
                                False,
                            )
                            getattr(self, to_add).update({gene.id: new_gene})
                    elif type(attr.values[0]) != str:
                        if gene.id in getattr(self, to_add).keys():
                            getattr(self, to_add).get(gene.id)._updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                gene.id,
                                gene.id,
                                model_id,
                                list(all_models_data.keys()),
                                False,
                            )
                            getattr(self, to_add).update({gene.id: new_gene})
                    else:
                        new_id = attr.values[0]
                        if new_id in self.assembly.keys():
                            self.assembly.get(new_id)._updateNewGene(gene.id, model_id)
                            getattr(self, model_id).update(
                                {new_id: self.assembly.get(new_id)}
                            )
                        else:
                            new_gene = NewGene(
                                new_id,
                                gene.id,
                                model_id,
                                list(all_models_data.keys()),
                                True,
                            )
                            self.assembly.update({new_id: new_gene})
                            getattr(self, model_id).update({new_id: new_gene})

    def __init__(self, all_models_data: dict, gene_folder, do_mix_conv_notconv: bool):
        self.assembly = {}
        for source in list(all_models_data.keys()):
            setattr(self, source, {})
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        if gene_folder is not None:
            self.__addNewGenes_conv(all_models_data, gene_folder, do_mix_conv_notconv)
        else:
            for model_id in list(all_models_data.keys()):
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id,
                            gene.id,
                            model_id,
                            list(all_models_data.keys()),
                            False,
                        )
                        self.assembly.update({gene.id: new_gene})
            for gene in self.assembly.values():
                for model_id in gene.in_models["models_list"]:
                    getattr(self, model_id).update({gene.id: gene})


class SuperModel:  # TODO REAL 30.08.23 add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """ Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values.  """

    def __find_connections(
        self,
        connection_knowledge: KnowledgeConnectingOldNew,
        all_models_data: dict,
        do_mix: bool,
        gene_folder,
    ):
        for met in self.metabolites.assembly.values():
            met._find_reactions(connection_knowledge)
        for r in self.reactions.assembly.values():
            r._find_reactants_products(connection_knowledge, "reactants")
            r._find_reactants_products(connection_knowledge, "products")
            r._find_metabolites(connection_knowledge)
            g_to_add = r._find_gene_and_gpr(connection_knowledge, gene_folder)
            for model_id, gene_ids in g_to_add.items():
                for g_id in gene_ids:
                    r.genes[model_id].append(self.genes.assembly[g_id])
        for gene in self.genes.assembly.values():
            gene._find_reactions(all_models_data, connection_knowledge)
        if not do_mix:
            for met in self.metabolites.notconverted.values():
                met._find_reactions(connection_knowledge, do_notconv=True)
            for r in self.reactions.notconverted.values():
                r._find_reactants_products(
                    connection_knowledge, "reactants", do_notconv=True
                )
                r._find_reactants_products(
                    connection_knowledge, "products", do_notconv=True
                )
                r._find_metabolites(connection_knowledge, do_notconv=True)
                g_to_add = r._find_gene_and_gpr(
                    connection_knowledge, gene_folder, do_notconv=True
                )
                for model_id, gene_ids in g_to_add.items():
                    for g_id in gene_ids:
                        if g_id in self.genes.notconverted.keys():
                            r.genes[model_id].append(self.genes.notconverted[g_id])
            for gene in self.genes.notconverted.values():
                gene._find_reactions(all_models_data, connection_knowledge, True)

    def __get_additional_attributes(
        self,
        model_ids: [str],
        connections: KnowledgeConnectingOldNew,
        do_mix_conv_notconv: bool,
    ):
        where_to_look = {"assembly": False}
        if not do_mix_conv_notconv:
            where_to_look.update({"notconverted": True})
        for atr, do_notconv in where_to_look.items():
            for met in getattr(self.metabolites, atr).values():
                for model_id in model_ids:
                    old_mets = connections.get_old_mets(model_id, met.id, do_notconv)
                    if old_mets:
                        met.formula.get(model_id).append(old_mets[0].formula)
            for r in getattr(self.reactions, atr).values():
                for mod_id in model_ids:
                    old_rs = connections.get_old_rs(mod_id, r.id, do_notconv)
                    if old_rs:
                        low_b = 0
                        upp_b = 0
                        subsys = []
                        for old_r in old_rs:
                            if old_r.lower_bound < low_b:
                                low_b = old_r.lower_bound
                            if old_r.upper_bound > upp_b:
                                upp_b = old_r.upper_bound
                            subsys.append(old_r.subsystem)
                        r.lower_bound.get(mod_id).append(low_b)
                        r.upper_bound.get(mod_id).append(upp_b)
                        r.subsystem.get(mod_id).append("#or#".join(subsys))

    def __swapReactantsAndProducts(self, r: NewReaction, sources_to_swap: list):
        for s in sources_to_swap:
            a = r.reactants.get(s)
            b = r.products.get(s)
            r.reactants[s] = b
            r.products[s] = a
            aa = r.lower_bound.get(s)[0] * -1
            bb = r.upper_bound.get(s)[0] * -1
            r.lower_bound[s] = [bb]
            r.upper_bound[s] = [aa]
            for met, koef in r.metabolites.get(s).items():
                r.metabolites.get(s)[met] = koef * -1

    def __runSwitchedMetabolites(self):
        for r in self.reactions.assembly.values():
            ex = False
            react_in = r.reactants[r.in_models["models_list"][0]]
            pro_in = r.products[r.in_models["models_list"][0]]
            for tmp in r.in_models["models_list"]:
                react_in = list(set(react_in) & set(r.reactants.get(tmp)))
                pro_in = list(set(pro_in) & set(r.products.get(tmp)))
                if (not r.reactants.get(tmp)) | (not r.products.get(tmp)):
                    ex = True
            if not ex:
                if (not react_in) | (not pro_in):
                    up = r.in_models["models_amount"] - 1
                    down = ceil(r.in_models["models_amount"] / 2) - 1
                    consist = []
                    for i in range(up, down, -1):
                        combinations = list(
                            itertools.combinations(r.in_models["models_list"], i)
                        )
                        for comb in combinations:
                            react_in_comb = r.reactants[comb[0]]
                            for c in comb:
                                react_in_comb = list(
                                    set(react_in_comb) & set(r.reactants.get(c))
                                )
                            if react_in_comb:
                                consist.append(comb)
                        if consist != []:
                            break
                    if len(consist) == 1:
                        # "Case 1: majority"
                        source_to_swap = list(
                            set(r.in_models["models_list"]) - set(consist[0])
                        )
                        self.__swapReactantsAndProducts(r, source_to_swap)
                    elif len(consist) == 2:
                        lb1 = 0
                        lb2 = 0
                        for tmp in r.in_models["models_list"]:
                            if tmp in consist[0]:
                                if r.lower_bound.get(tmp)[0] < lb1:
                                    lb1 = r.lower_bound.get(tmp)[0]
                            if tmp in consist[1]:
                                if r.lower_bound.get(tmp)[0] < lb2:
                                    lb2 = r.lower_bound.get(tmp)[0]
                        swap = None
                        if (lb1 >= 0) & (lb2 < 0):
                            swap = consist[1]
                        if (lb1 < 0) & (lb2 >= 0):
                            swap = consist[0]
                        if swap:
                            # "Case 2: boundary"
                            self.__swapReactantsAndProducts(r, swap)
                        else:
                            # "Case 3: Nothing sort"
                            sel = sorted(r.in_models["models_list"])[0]
                            not_sel = []
                            for tmp in sorted(r.in_models["models_list"])[1:]:
                                if not (
                                    set(r.reactants.get(tmp))
                                    & set(r.reactants.get(sel))
                                ):
                                    not_sel.append(tmp)
                            self.__swapReactantsAndProducts(r, not_sel)
                    # len(consist) is expected to be only 1 or 2
                    else:
                        warnings.warn(
                            f"Warning! Something went wrong with swaping "
                            f"metabolites for {r.id}."
                            f"Can enter consist more 2."
                            f"Reactants: {r.reactants}. "
                            f"Products: {r.products}. Consist: {consist}"
                        )
                        # "Case 3: Nothing sort"
                        sel = sorted(r.in_models["models_list"])[0]
                        not_sel = []
                        for tmp in sorted(r.in_models["models_list"])[1:]:
                            if not (
                                set(r.reactants.get(tmp)) & set(r.reactants.get(sel))
                            ):
                                not_sel.append(tmp)
                        self.__swapReactantsAndProducts(r, not_sel)

    def __assemble_attributes(self, and_as_solid: bool, do_mix_conv_notconv: bool):
        where_assemble = ["assembly"]
        if not do_mix_conv_notconv:
            where_assemble.append("notconverted")
        for atr in where_assemble:
            for met in getattr(self.metabolites, atr).values():
                ass_r = getCoreConnections(met.reactions, 1, operator.ge, self.sources)
                met.reactions.update({"assembly": ass_r})
            for gene in getattr(self.genes, atr).values():
                ass_rg = getCoreConnections(
                    gene.reactions, 1, operator.ge, self.sources
                )
                gene.reactions.update({"assembly": ass_rg})
            for react in getattr(self.reactions, atr).values():
                ass_reactants = getCoreConnections(
                    react.reactants, 1, operator.ge, self.sources
                )
                ass_products = getCoreConnections(
                    react.products, 1, operator.ge, self.sources
                )
                ass_genes = getCoreConnections(
                    react.genes, 1, operator.ge, self.sources
                )
                ass_gpr = getCoreGPR(
                    react.gene_reaction_rule,
                    1,
                    operator.ge,
                    self.sources,
                    and_as_solid,
                )
                ass_lower_bound = getCoreLowerBounds(
                    react.lower_bound, 1, react.in_models["models_list"]
                )
                ass_upper_bound = getCoreUpperBounds(
                    react.upper_bound, 1, react.in_models["models_list"]
                )
                react.reactants.update({"assembly": ass_reactants})
                react.products.update({"assembly": ass_products})
                react.genes.update({"assembly": ass_genes})
                react.gene_reaction_rule.update({"assembly": ass_gpr})
                react.lower_bound.update({"assembly": ass_lower_bound})
                react.upper_bound.update({"assembly": ass_upper_bound})
                core_metabolites = getCoreCoefficients(
                    react.metabolites,
                    react.reactants,
                    react.products,
                    "assembly",
                    1,
                    react.in_models["models_list"],
                )
                react.metabolites.update({"assembly": core_metabolites})

    def __init__(
        self,
        final_m_sel: dict,
        final_m_not_sel: dict,
        final_r_sel: dict,
        final_r_not_sel: dict,
        all_models_data: dict,
        additional_periplasmic_m: dict,
        periplasmic_r: dict,
        m_db_info: pd.core.frame.DataFrame,
        r_db_info: pd.core.frame.DataFrame,
        gene_folder,
        do_mix_conv_notconv: bool,
        and_as_solid: bool,
    ):
        self.sources = list(all_models_data.keys())
        self.metabolites = SetofNewElements(
            "metabolites",
            final_m_sel,
            final_m_not_sel,
            self.sources,
            m_db_info,
            do_mix_conv_notconv,
            additional_periplasmic_m,
        )
        self.reactions = SetofNewElements(
            "reactions",
            final_r_sel,
            final_r_not_sel,
            self.sources,
            r_db_info,
            do_mix_conv_notconv,
        )
        self.genes = SetofNewGenes(all_models_data, gene_folder, do_mix_conv_notconv)

        if do_mix_conv_notconv:
            final_m_all = defaultdict(dict)
            for model_id in final_m_sel.keys():
                final_m_all[model_id] = (
                    final_m_sel[model_id] | final_m_not_sel[model_id]
                )
            m_connection_dicts = self.metabolites._makeForwardBackward(
                all_models_data,
                final_m_all,
                "metabolites",
                additional=additional_periplasmic_m,
            )
            final_r_all = defaultdict(dict)
            for model_id in final_r_sel.keys():
                final_r_all[model_id] = (
                    final_r_sel[model_id] | final_r_not_sel[model_id]
                )
            r_connection_dicts = self.reactions._makeForwardBackward(
                all_models_data, final_r_all, "reactions",
            )
        else:
            m_connection_dicts = self.metabolites._makeForwardBackward(
                all_models_data,
                final_m_sel,
                "metabolites",
                additional=additional_periplasmic_m,
                not_selected=final_m_not_sel,
            )
            r_connection_dicts = self.reactions._makeForwardBackward(
                all_models_data, final_r_sel, "reactions", not_selected=final_r_not_sel
            )
        connection_knowledge = KnowledgeConnectingOldNew(
            m_connection_dicts,
            r_connection_dicts,
            additional_periplasmic_m,
            periplasmic_r,
            gene_folder,
        )
        self.__find_connections(
            connection_knowledge, all_models_data, do_mix_conv_notconv, gene_folder,
        )
        self.__get_additional_attributes(
            self.sources, connection_knowledge, do_mix_conv_notconv
        )
        self.__runSwitchedMetabolites()
        self.__assemble_attributes(and_as_solid, do_mix_conv_notconv)

    def get_short_name_len(self) -> int:
        for i in range(len(max(self.sources, key=len)) + 1):
            short = []
            for source in self.sources:
                short.append(source[:i])
            if len(set(short)) == len(self.sources):
                return i

    def at_least_in(self, number_of_model: int, and_as_solid=False):
        if (
            type(number_of_model) != int
            or number_of_model < 1
            or number_of_model > len(self.sources)
        ):
            raise ValueError("Number to check does not fit the number of models")
        elif number_of_model == 1:
            raise ValueError(
                "Features in at least 1 model are already found in assembly. "
                "You do not need to run this comparison separately"
            )
        else:
            getCore(self, number_of_model, operator.ge, and_as_solid)

    def exactly_in(self, number_of_model: int, and_as_solid=False):
        if (
            type(number_of_model) != int
            or number_of_model < 1
            or number_of_model > len(self.sources)
        ):
            raise ValueError("Number to check does not fit the number of models")
        else:
            getCore(self, number_of_model, operator.eq, and_as_solid)

    def present(self, yes=None, no=None, short_name_len=None, and_as_solid=False):
        if yes is None and no is None:
            raise ValueError(
                "Both models present and models not present are not provided. "
                "Please provide at least one of the list"
            )
        elif (yes is not None and type(yes) != list) or (
            no is not None and type(no) != list
        ):
            raise ValueError(
                "Present or not present models are in wrong type. "
                "Please provide lists"
            )
        else:
            if yes is None:
                yes = []
            if no is None:
                no = []
            wrong_yes = set(yes) - set(self.sources)
            wrong_no = set(no) - set(self.sources)
            if wrong_yes or wrong_no:
                raise ValueError(
                    f"Some of input models are not in supermodel. "
                    f"Maybe {wrong_yes} or {wrong_no}"
                    f"Please check the input ids"
                )
            else:
                if short_name_len is None:
                    short_name_len = self.get_short_name_len()
                getDifference(self, yes, no, and_as_solid, short_name_len)

    def get_venn_segments(self, short_name_len=None, and_as_solid=False):
        """ Getting metabolites and reactions networks for each Venn segment in Venn
        diagram."""
        if short_name_len is None:
            short_name_len = self.get_short_name_len()
        combinations = []
        for i in range(1, len(self.sources)):
            combinations.extend(itertools.combinations(self.sources, i))
        for combo in combinations:
            yes = sorted(list(combo))
            no = sorted((list(set(self.sources) - set(combo))))
            getDifference(self, yes, no, and_as_solid, short_name_len)

    def get_intersection(self, and_as_solid=False):
        getCore(self, len(self.sources), operator.ge, and_as_solid)

    def get_all_confident_levels(self, and_as_solid=False):
        for i in range(len(self.sources), 1, -1):
            self.at_least_in(i, and_as_solid=and_as_solid)

    def write_supermodel_to_pkl(self, output_name: str, recursion_limit=None):
        if not output_name.endswith(".pkl"):
            raise ValueError("Wrong extension of the file")
        if exists(output_name):
            raise ValueError("File already exist, change the name")
        else:
            max_rec = 0x100000
            resource.setrlimit(
                resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY]
            )
            sys.setrecursionlimit(max_rec)
            with open(output_name, "wb") as fh:
                dill.dump(self, fh)


def read_supermodel_from_pkl(input_name: str):
    supermodel = dill.load(open(input_name, "rb"))
    return supermodel
