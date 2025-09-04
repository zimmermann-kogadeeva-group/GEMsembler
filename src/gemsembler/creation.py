import itertools

# import resource #Windows_Fix
import json
import operator
import re
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
        do_old_genes: dict,
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
        self.do_genes = {}
        for model_id in self.m_go_old_new.keys():
            if gene_folder is None:
                self.g_conversion_tables[model_id] = None
                self.do_genes[model_id] = do_old_genes[model_id]
                continue
            blast_file = gene_folder / (model_id + "_blast.tsv")
            try:
                conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            except:
                self.g_conversion_tables[model_id] = None
                self.do_genes[model_id] = do_old_genes[model_id]
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
                self.do_genes[model_id] = True

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
    """New object class - one metabolite or reaction for supermodel."""

    def __init__(
        self,
        is_loading,
        args_dict,  # structure: {"Type": "NameOfClass", "args":{"new_ids":......}}
        #       new_id: str, # Old input
        #       old_id: str,#
        #       compartments: [str],
        #       source: str,
        #       possible_sources: [str],
        #       converted: bool,
    ):
        if args_dict["type"] != "NewElement":  # Houston, we have a problem
            raise ValueError(
                f"The args for the NewElement objects are from different object ({args_dict['type']})!"
            )
            # warnings.warn(f"") # for warning

        args = args_dict["args"]

        if is_loading:
            self.id = args["new_id"]
            self.compartments = args["compartments"]
            self.sources = args["sources"]
            self.in_models = args["in_models"]
            self.annotation = args["annotation"]
            self.converted = args["converted"]
            return

        # else: initial init
        self.id = args["new_id"]
        self.compartments = {"assembly": args["compartments"]}
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [args["source"]]}
        self.annotation = {}
        self.converted = args["converted"]
        for ps in args["possible_sources"]:
            if ps == args["source"]:
                self.compartments.update({ps: args["compartments"]})
                self.sources.update({ps: 1})
                self.annotation.update({ps: [args["old_id"]]})
            else:
                self.compartments.update({ps: []})
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})
        return

    def _args_to_dict(self):  # GGE
        out_json_dict = {}
        out_json_dict["new_id"] = self.id
        out_json_dict["compartments"] = self.compartments
        out_json_dict["sources"] = self.sources
        out_json_dict["in_models"] = self.in_models
        out_json_dict["annotation"] = self.annotation
        out_json_dict["converted"] = self.converted
        return out_json_dict

    def _to_json(self):  # GGE
        return json.dumps(self._args_to_dict())

    def _get_replace_tag(self):  # GGE
        tag = ""
        if self.converted:
            tag = "~"
        else:
            tag = "%"  # unique tag for not converted, as they could have similar id as converted one
        return f"{tag}{self.id}"

    def _update_new_element(
        self,
        id_to_update: str,
        compart_to_update: [str],
        source: str,
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
        is_loading,
        args_dict,
        # structure: {"type": "NewMetaboliteOrReaction", "args":{"new_ids": , "old_id": , "compartments": , "source": , "possible_sources", "converted": , "db_info": }}
        #        new_id: str, # Old input
        #        old_id: str,
        #        compartments: [str],
        #        source: str,
        #        possible_sources: [str],
        #        converted: bool,
        #        m_database_info: pd.core.frame.DataFrame, -> args["db_info"]
    ):
        if (
            args_dict["type"] != "NewMetaboliteOrReaction"
            and args_dict["type"] != "NewMetabolite"
        ):  # Houston, we have a problem
            raise ValueError(
                f"The args for the NewMetabolite objects are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        super().__init__(is_loading, {"type": "NewElement", "args": args})

        if is_loading:
            self.name = args["name"]
            self.reactions = args["reactions"]
            self.formula = args["formula"]
            if "formula_bigg" in args.keys():
                self.formula_bigg = args["formula_bigg"]
            else:
                self.formula_bigg = None
            if "charge" in args.keys():
                self.charge = args["charge"]
            else:
                self.charge = None
            if "charge_bigg" in args.keys():
                self.charge_bigg = args["charge_bigg"]
            else:
                self.charge_bigg = None
            return

        # else: initial init
        if args["converted"]:
            all_comp_bigg = (
                args["db_info"].bigg_id.str.slice(-1).drop_duplicates().pipe("".join)
            )

            # Get specific row from db_info table
            id_noc = re.sub(f"_([{all_comp_bigg}])$", "", args["new_id"])
            row = args["db_info"].query(f"universal_bigg_id == '{id_noc}'")

            if row.empty:
                raise RuntimeError(
                    f"Issue with `{id_noc}` "
                    "(ID with compartment: {args['new_id']}) "
                    "- not found in the table."
                )

            # Get additional attributes from db_info table if present
            name = row.name.iloc[0] if "name" in row else None
            formula = row.formula.iloc[0] if "formula" in row else None
            charge = row.charge.iloc[0] if "charge" in row else None
        else:
            name = "Not converted"
            formula = None
            charge = None
        self.name = name
        self.formula_bigg = formula
        self.charge_bigg = charge
        self.formula = {k: [] for k in args["possible_sources"]}
        self.charge = {k: [] for k in args["possible_sources"]}
        self.reactions = {k: [] for k in args["possible_sources"]}
        self.reactions.update({"assembly": [], "comparison": {}})

    def _args_to_dict(self):  # GGE
        out_json_dict = super()._args_to_dict()
        out_json_dict["name"] = self.name
        out_json_dict["formula"] = self.formula
        out_json_dict["charge"] = self.charge
        out_json_dict["formula_bigg"] = self.formula_bigg
        out_json_dict["charge_bigg"] = self.charge_bigg
        out_json_dict["reactions"] = self.reactions
        out_json_dict["type"] = "NewMetabolite"
        return out_json_dict

    def _to_json(self):  # GGE
        args_dict_out = self._args_to_dict()
        json_dict_out = {}
        keys_with_problems = ["reactions"]
        for key in args_dict_out.keys():
            if key not in keys_with_problems:
                json_dict_out[key] = args_dict_out[key]

        # Rewriting pointers to objects with just their names (ids?)
        json_dict_out["reactions"] = {}

        for c_key, c_val in args_dict_out["reactions"].items():
            if not c_val:  # dict or list is empty
                json_dict_out["reactions"][c_key] = c_val
                continue
            elif type(c_val) is list:
                json_dict_out["reactions"][c_key] = []
                for l in c_val:
                    json_dict_out["reactions"][c_key].append(f"{l._get_replace_tag()}")
                continue
            elif c_key == "comparison":
                json_dict_out["reactions"][c_key] = {}
            elif type(c_val) is dict:
                raise Warning(f"Should not be dict {c_val} for {c_key}")
                # json_dict_out["reactions"][c_key] = {}
                # for c_key2, c_val2 in c_val.items():
                #     json_dict_out["reactions"][c_key][c_key2] = []
                #     for el in c_val2:
                #         json_dict_out["reactions"][c_key][c_key2].append(
                #             f"{el._get_replace_tag()}"
                #         )

        return json.dumps(json_dict_out)

    def _replace_tags_with_objects(self, all_objects_by_id: dict):  # GGE
        new_reactions = {}
        for c_key, c_val in self.reactions.items():
            if not c_val:  # dict or list is empty
                new_reactions[c_key] = c_val
            elif type(c_val) is list:
                new_reactions[c_key] = []
                for l in range(len(c_val)):
                    if c_val[l] in all_objects_by_id:
                        new_reactions[c_key].append(all_objects_by_id[c_val[l]])
                    else:
                        new_reactions[c_key].append(c_val[l])
            elif type(c_val) is dict:
                new_reactions[c_key] = {}
                for c_key2, c_val2 in c_val.items():
                    new_reactions[c_key][c_key2] = []
                    for el in range(len(c_val2)):
                        if c_val2[el] in all_objects_by_id:
                            new_reactions[c_key][c_key2].append(
                                all_objects_by_id[c_val2[el]]
                            )
                        else:
                            new_reactions[c_key][c_key2].append(c_val2[el])
        self.reactions = new_reactions

        return

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
        is_loading,
        args_dict,
        # structure: {"type": "NewMetaboliteOrReaction", "args":{"new_ids": , "old_id": , "compartments": , "source": , "possible_sources", "converted": , "db_info": }}
        #        new_id: str,  # Old input
        #        old_id: str,
        #        compartments: [str],
        #        source: str,
        #        possible_sources: [str],
        #        converted: bool,
        #        r_database_info: pd.core.frame.DataFrame,  -> args["db_info"]
    ):
        if (
            args_dict["type"] != "NewMetaboliteOrReaction"
            and args_dict["type"] != "NewReaction"
        ):  # Houston, we have a problem
            raise ValueError(
                f"The args for the NewReaction objects are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        super().__init__(  # GGE template for testing - continue later
            is_loading, {"type": "NewElement", "args": args}
        )

        if is_loading:
            self.name = args["name"]
            self.reaction = args["reaction"]
            self.reactants = args["reactants"]
            self.products = args["products"]
            self.metabolites = args["metabolites"]
            self.lower_bound = args["lower_bound"]
            self.upper_bound = args["upper_bound"]
            self.subsystem = args["subsystem"]
            self.genes = args["genes"]
            self.gene_reaction_rule = args["gene_reaction_rule"]
            return

        # else: initial init
        if args["converted"]:
            id_noc = args["new_id"].replace("sink_", "DM_")
            name = args["db_info"][args["db_info"]["bigg_id"] == id_noc]["name"]
            if (not name.empty) and (not name.isnull().values.any()):
                name = name.values[0]
            else:
                name = ""
            equation = args["db_info"][args["db_info"]["bigg_id"] == id_noc][
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
        base_keys = args["possible_sources"] + ["assembly"]
        self.reactants = {k: [] for k in base_keys}
        self.reactants.update({"comparison": {}})
        self.products = {k: [] for k in base_keys}
        self.products.update({"comparison": {}})
        self.metabolites = {k: {} for k in base_keys + ["comparison"]}
        self.lower_bound = {k: [] for k in base_keys}
        self.lower_bound.update({"comparison": {}})
        self.upper_bound = {k: [] for k in base_keys}
        self.upper_bound.update({"comparison": {}})
        self.subsystem = {k: [] for k in args["possible_sources"]}
        self.genes = {k: [] for k in base_keys}
        self.genes.update({"comparison": {}})
        self.gene_reaction_rule = {k: [] for k in base_keys}
        self.gene_reaction_rule.update(
            {k + "_mixed": [] for k in args["possible_sources"]}
        )
        self.gene_reaction_rule.update({"comparison": {}})

    def _args_to_dict(self):  # GGE
        out_json_dict = super()._args_to_dict()
        out_json_dict["name"] = self.name
        out_json_dict["reaction"] = self.reaction
        out_json_dict["lower_bound"] = {**self.lower_bound, "comparison": {}}
        out_json_dict["upper_bound"] = {**self.upper_bound, "comparison": {}}
        out_json_dict["subsystem"] = self.subsystem
        out_json_dict["gene_reaction_rule"] = {
            **self.gene_reaction_rule,
            "comparison": {},
        }
        # DICTS WITH LINKS TO OBJECTS
        out_json_dict["metabolites"] = self.metabolites
        out_json_dict["products"] = self.products
        out_json_dict["reactants"] = self.reactants
        out_json_dict["genes"] = self.genes
        out_json_dict["type"] = "NewReaction"
        return out_json_dict

    def _to_json(self):  # GGE
        # Can't just say json_dict_out = args_dict_out and than edit json_dict_out - we will also edit the args_dict_out dict
        # To solve this pointer problem, I chose the folowing way (GGE)
        args_dict_out = self._args_to_dict()
        json_dict_out = {}
        keys_with_problems = ["reactants", "metabolites", "genes", "products"]
        for key in args_dict_out.keys():
            if key not in keys_with_problems:
                json_dict_out[key] = args_dict_out[key]

        # Rewriting pointers to objects with just their ids
        json_dict_out["reactants"] = {}
        json_dict_out["metabolites"] = {}
        json_dict_out["genes"] = {}
        json_dict_out["products"] = {}
        for a_key, a_val in args_dict_out[
            "reactants"
        ].items():  # vals are dicts or lists
            if not a_val:  # dict or list is empty
                json_dict_out["reactants"][a_key] = a_val
            elif type(a_val) is list:
                json_dict_out["reactants"][a_key] = []
                for elm in a_val:
                    json_dict_out["reactants"][a_key].append(
                        f"{elm._get_replace_tag()}"
                    )
            elif a_key == "comparison":
                json_dict_out["reactants"][a_key] = {}
            elif type(a_val) is dict:
                raise Warning(f"Should not be dict {a_val} for {a_key}")
                # json_dict_out["reactants"][a_key] = {}
                # for child_key in a_val.keys():
                #     json_dict_out["reactants"][a_key][child_key] = []
                #     for elem in a_val[child_key]:
                #         json_dict_out["reactants"][a_key][child_key].append(
                #             f"{elem._get_replace_tag()}"
                #         )
            else:
                raise ValueError(
                    f"NewReaction.to_json: Unexpected type of value for 'reactants'! ({type(a_val)})!"
                )

        for a_key, a_val in args_dict_out[
            "metabolites"
        ].items():  # vals are dicts or lists
            if not a_val:  # dict or list is empty
                json_dict_out["metabolites"][a_key] = a_val
            elif type(a_val) is dict:
                json_dict_out["metabolites"][a_key] = {}
                if a_key == "comparison":
                    continue
                lockal_dict = json_dict_out["metabolites"][a_key]
                for (
                    a_key2,
                    a_val2,
                ) in a_val.items():  # GGE here the key is object, and value is normal
                    new_key2 = ""
                    if type(a_key2) is not NewMetabolite:
                        new_key2 = a_key2
                    else:
                        new_key2 = f"{a_key2._get_replace_tag()}"
                    if type(a_val2) is not dict:
                        lockal_dict[new_key2] = a_val2
                    else:
                        raise Warning("No comparison currently")
                        # lockal_dict[new_key2] = {}
                        # for a_key3, a_val3 in a_val2.items():
                        #     if type(a_key3) is not NewMetabolite:
                        #         raise ValueError(
                        #             f"NewReaction.to_json: Unexpected type of value for 'metabolites'! a_key3 ({type(a_key3)})!"
                        #         )
                        #     new_key3 = f"{a_key3._get_replace_tag()}"
                        #     lockal_dict[new_key2][new_key3] = a_val3
            else:
                raise ValueError(
                    f"NewReaction.to_json: Unexpected type of value for 'metabolites'! ({type(a_val)})!"
                )

        for c_key, c_val in args_dict_out["genes"].items():
            if not c_val:  # dict or list is empty
                json_dict_out["genes"][c_key] = c_val
                continue
            elif type(c_val) is list:
                json_dict_out["genes"][c_key] = []
                for l in c_val:
                    json_dict_out["genes"][c_key].append(f"{l._get_replace_tag()}")
                continue
            elif c_key == "comparison":
                json_dict_out["genes"][c_key] = {}
            elif type(c_val) is dict:
                raise Warning(f"Should not b dict {c_val} for {c_key}")
                # json_dict_out["genes"][c_key] = {}
                # for c_key2, c_val2 in c_val.items():
                #     json_dict_out["genes"][c_key][c_key2] = []
                #     #                    print(f"c_key2 is {c_key2} c_val2 is {c_val2}")
                #     for el in c_val2:
                #         json_dict_out["genes"][c_key][c_key2].append(
                #             f"{el._get_replace_tag()}"
                #         )

        for c_key, c_val in args_dict_out["products"].items():
            if not c_val:  # dict or list is empty
                json_dict_out["products"][c_key] = c_val
                continue
            elif type(c_val) is list:
                json_dict_out["products"][c_key] = []
                for l in c_val:
                    json_dict_out["products"][c_key].append(f"{l._get_replace_tag()}")
                continue
            elif c_key == "comparison":
                json_dict_out["genes"][c_key] = {}
            elif type(c_val) is dict:
                raise Warning(f"Should not b dict {c_val} for {c_key}")
                # json_dict_out["products"][c_key] = {}
                # for c_key2, c_val2 in c_val.items():
                #     json_dict_out["products"][c_key][c_key2] = []
                #     #                    print(f"c_key2 is {c_key2} c_val2 is {c_val2}")
                #     for el in c_val2:
                #         json_dict_out["products"][c_key][c_key2].append(
                #             f"{el._get_replace_tag()}"
                #         )

        return json.dumps(json_dict_out)

    def _replace_tags_with_objects(self, all_objects_by_id: dict):  # GGE
        for a_key, a_val in self.reactants.items():  # vals are dicts or lists
            if not a_val:  # dict or list is empty
                continue
            elif type(a_val) is list:
                for elm in range(len(a_val)):
                    if a_val[elm] in all_objects_by_id.keys():
                        a_val[elm] = all_objects_by_id[a_val[elm]]
            elif type(a_val) is dict:
                for child_key in a_val.keys():
                    for elem in range(len(a_val[child_key])):
                        if a_val[child_key][elem] in all_objects_by_id.keys():
                            a_val[child_key][elem] = all_objects_by_id[
                                a_val[child_key][elem]
                            ]

        new_metabolites = {}
        for a_key, a_val in self.metabolites.items():  # vals are dicts or lists
            if not a_val:  # dict or list is empty
                new_metabolites[a_key] = a_val
                continue
            elif type(a_val) is dict:
                new_metabolites[a_key] = {}
                for (
                    a_key2,
                    a_val2,
                ) in a_val.items():  # GGE here the key is object, and value is normal
                    if type(a_val2) is not dict:
                        if a_key2 in all_objects_by_id.keys():
                            new_metabolites[a_key][all_objects_by_id[a_key2]] = a_val[
                                a_key2
                            ]
                        else:
                            new_metabolites[a_key][a_key2] = a_val[a_key2]
                    else:
                        new_metabolites[a_key][a_key2] = {}
                        for a_key3, a_val3 in a_val2.items():
                            if a_key3 in all_objects_by_id.keys():
                                new_metabolites[a_key][a_key2][
                                    all_objects_by_id[a_key3]
                                ] = a_val2[a_key3]
                            else:
                                new_metabolites[a_key][a_key2][a_key3] = a_val2[a_key3]
        self.metabolites = new_metabolites

        for c_key, c_val in self.genes.items():
            if not c_val:  # dict or list is empty
                continue
            elif type(c_val) is list:
                for l in range(len(c_val)):
                    if c_val[l] in all_objects_by_id.keys():
                        c_val[l] = all_objects_by_id[c_val[l]]
            elif type(c_val) is dict:
                for c_key2, c_val2 in c_val.items():
                    for e in range(len(c_val2)):
                        if c_val2[e] in all_objects_by_id.keys():
                            c_val2[e] = all_objects_by_id[c_val2[e]]

        for c_key, c_val in self.products.items():
            if not c_val:  # dict or list is empty
                continue
            elif type(c_val) is list:
                for l in range(len(c_val)):
                    if c_val[l] in all_objects_by_id.keys():
                        c_val[l] = all_objects_by_id[c_val[l]]
            elif type(c_val) is dict:
                for c_key2, c_val2 in c_val.items():
                    for e in range(len(c_val2)):
                        if c_val2[e] in all_objects_by_id.keys():
                            c_val2[e] = all_objects_by_id[c_val2[e]]
        return

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
        self,
        connections: KnowledgeConnectingOldNew,
        m_type: str,
        do_notconv=False,
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
                    old_react_react_prod,
                    model_id,
                    connections,
                    do_notconv,
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
        self,
        connections: KnowledgeConnectingOldNew,
        do_notconv=False,
    ):
        genes_to_add = defaultdict(list)
        for model_id in self.in_models["models_list"]:
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
                    if connections.do_genes[model_id]:
                        genes_to_add[model_id].append(pot_new_g_id)
                old_gpr = oldr.gene_reaction_rule
                new_gpr, mix_gpr = makeNewGPR(old_gpr, gene_convert)
                if new_gpr:
                    new_gpr_unite_r.append(new_gpr)
                self.gene_reaction_rule.get(model_id + "_mixed").append(mix_gpr)
            if len(new_gpr_unite_r) == 1 and connections.do_genes[model_id]:
                self.gene_reaction_rule.get(model_id).append(new_gpr_unite_r[0])
            elif len(new_gpr_unite_r) >= 1 and connections.do_genes[model_id]:
                united_gpr = uniteGPR(new_gpr_unite_r)
                self.gene_reaction_rule.get(model_id).append(united_gpr)
        return genes_to_add


class SetofNewElements:
    """Setting dictionaries for all metabolites or reactions:
    selected for supermodel - self.assembly and not selected - self.notconverted."""

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
                                False,
                                {
                                    "type": "NewMetaboliteOrReaction",
                                    "args": {
                                        "new_id": new_id + "_convert_" + str(convered),
                                        "old_id": key,
                                        "compartments": comp,
                                        "source": mod_id,
                                        "possible_sources": model_ids,
                                        "converted": convered,
                                        "db_info": db_info,
                                    },
                                },
                            )
                            getattr(self, where_to_add).update(
                                {new_id + "_convert_" + str(convered): new_conv}
                            )
                    else:
                        new = new_elements[element_type](
                            False,
                            {
                                "type": "NewMetaboliteOrReaction",
                                "args": {
                                    "new_id": new_id,
                                    "old_id": key,
                                    "compartments": comp,
                                    "source": mod_id,
                                    "possible_sources": model_ids,
                                    "converted": convered,
                                    "db_info": db_info,
                                },
                            },
                        )
                        getattr(self, where_to_add).update({new_id: new})

    def __init__(
        self,
        is_loading: bool,
        args_dict,
        # structure: {"type": "SetofNewElements", "args": {"element_type": , "selected": ,  "not_selected": ,  "model_ids": ,  "db_info": ,  "do_mix_conv_notconv": }}
        additional=None,
        #        element_type: str, # Old input
        #        selected: dict,
        #        not_selected: dict,
        #        model_ids: [str],
        #        db_info: pd.core.frame.DataFrame,
        #        do_mix_conv_notconv: bool,
    ):
        if args_dict["type"] != "SetofNewElements":  # Houston, we have a problem
            raise ValueError(
                f"The args for the SetofNewElements objects are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        if is_loading:
            print("loading SetofNewElements...")
            args = json.loads(args)
            self.assembly = {}  # args["assembly"]
            self.comparison = defaultdict(dict)  # args["comparison"]
            self.notconverted = {}  # args["notconverted"]

            for a_key, a_val in args["assembly"].items():
                object_dict = json.loads(a_val)
                if object_dict["type"] == "NewReaction":
                    self.assembly.update(
                        {
                            a_key: NewReaction(
                                True, {"type": "NewReaction", "args": object_dict}
                            )
                        }
                    )
                elif object_dict["type"] == "NewMetabolite":
                    self.assembly.update(
                        {
                            a_key: NewMetabolite(
                                True, {"type": "NewMetabolite", "args": object_dict}
                            )
                        }
                    )

            # for a_key, a_val in args["comparison"].items():
            #     for a_key2, a_val2 in a_val.items():
            #         object_dict = json.loads(a_val2)
            #         if object_dict["type"] == "NewReaction":
            #             self.comparison[a_key].update(
            #                 {
            #                     a_key2: NewReaction(
            #                         True, {"type": "NewReaction", "args": object_dict}
            #                     )
            #                 }
            #             )
            #         elif object_dict["type"] == "NewMetabolite":
            #             self.comparison[a_key].update(
            #                 {
            #                     a_key2: NewMetabolite(
            #                         True, {"type": "NewMetabolite", "args": object_dict}
            #                     )
            #                 }
            #             )

            for a_key, a_val in args["notconverted"].items():
                object_dict = json.loads(a_val)
                if object_dict["type"] == "NewReaction":
                    self.notconverted.update(
                        {
                            a_key: NewReaction(
                                True, {"type": "NewReaction", "args": object_dict}
                            )
                        }
                    )
                elif object_dict["type"] == "NewMetabolite":
                    self.notconverted.update(
                        {
                            a_key: NewMetabolite(
                                True, {"type": "NewMetabolite", "args": object_dict}
                            )
                        }
                    )
            return

        self.assembly = {}
        for source in args["selected"].keys():
            setattr(self, source, {})
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        self.__add_new_elements(
            args["element_type"],
            args["selected"],
            "assembly",
            args["model_ids"],
            True,
            args["db_info"],
        )
        if additional:
            self.__add_new_elements(
                args["element_type"],
                additional,
                "assembly",
                args["model_ids"],
                True,
                args["db_info"],
            )
        if args["do_mix_conv_notconv"]:
            self.__add_new_elements(
                args["element_type"],
                args["not_selected"],
                "assembly",
                args["model_ids"],
                False,
                args["db_info"],
            )
        else:
            self.__add_new_elements(
                args["element_type"],
                args["not_selected"],
                "notconverted",
                args["model_ids"],
                False,
                args["db_info"],
            )
        for new_id, new_obj in self.assembly.items():
            for model_id in new_obj.in_models["models_list"]:
                getattr(self, model_id).update({new_id: new_obj})

    def _args_to_dict(self):  # GGE
        out_dict = {}
        out_dict["assembly"] = self.assembly
        out_dict["comparison"] = {}  # self.comparison
        out_dict["notconverted"] = self.notconverted
        return out_dict

    def _to_json(self):  # GGE
        args_dict_out = self._args_to_dict()
        json_dict_out = {"assembly": {}, "comparison": {}, "notconverted": {}}
        # rewrite each element as json, than no problem writing all args
        for a_key, a_val in args_dict_out["assembly"].items():
            json_dict_out["assembly"].update({a_key: a_val._to_json()})
        # for a_key, a_val in args_dict_out["comparison"].items():
        #     for a_key2, a_val2 in a_val.items():
        #         a_val.update({a_key2: a_val2._to_json()})
        for a_key, a_val in args_dict_out["notconverted"].items():
            json_dict_out["notconverted"].update({a_key: a_val._to_json()})
        return json.dumps(json_dict_out)

    def _generate_objects_dict(self):  # GGE
        all_objects_by_id = {}

        for a_key, a_val in self.assembly.items():
            if a_val._get_replace_tag() in all_objects_by_id:
                print("Already exists!!!! NewEl, assembly")
            all_objects_by_id[a_val._get_replace_tag()] = a_val

        #        comparison - can be ignored (no unique elements)

        for a_key, a_val in self.notconverted.items():
            if a_val._get_replace_tag() in all_objects_by_id:
                print(f"Already exists!!!! NewEl, notconverted: {a_val.id}")
            all_objects_by_id[a_val._get_replace_tag()] = a_val

        return all_objects_by_id

    def _replace_tilda_tags(self, all_objects_by_id: dict):  # GGE

        for a_key, a_val in self.assembly.items():
            a_val._replace_tags_with_objects(all_objects_by_id)

        for a_key, a_val in self.comparison.items():
            for a_key2, a_val2 in a_val.items():
                a_val2._replace_tags_with_objects(all_objects_by_id)

        for a_key, a_val in self.notconverted.items():
            a_val._replace_tags_with_objects(all_objects_by_id)

        return

    def _makeForwardBackward(
        self,
        all_models: dict,
        selected: dict,
        obj_type: "metabolites" or "reactions",
        additional=None,
        not_selected=None,
    ):
        """Creating dictionaries linking metabolites/reactions:
        NewObject in supermodel with old original ID and OldObject in original models with new ID in supermodel
        """
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
        is_loading,
        args_dict,
        # structure: {"type": "NewGene", "args":{"new_id": , "old_id": , "source": , "possible_sources", "converted": }}
        #        new_id: str, # Old input
        #        old_id: str,
        #        source: str,
        #        possible_sources: [str],
        #        converted: bool,
    ):
        if args_dict["type"] != "NewGene":  # Houston, we have a problem
            raise ValueError(
                f"The args for the NewGene objects are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        if is_loading:
            self.id = args["new_id"]
            self.sources = args["sources"]
            self.converted = args["converted"]
            self.in_models = args["in_models"]
            self.annotation = args["annotation"]
            self.reactions = args["reactions"]
            return

        # else: initial init
        self.id = args["new_id"]
        self.sources = {}
        self.converted = args["converted"]
        self.in_models = {"models_amount": 1, "models_list": [args["source"]]}
        self.annotation = {}
        self.reactions = {"assembly": [], "comparison": {}}
        for ps in args["possible_sources"]:
            self.reactions.update({ps: []})
            if ps == args["source"]:
                self.sources.update({ps: 1})
                self.annotation.update({ps: [args["old_id"]]})
            else:
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def _args_to_dict(self):  # GGE
        out_json_dict = {}
        out_json_dict["new_id"] = self.id
        out_json_dict["sources"] = self.sources
        out_json_dict["converted"] = self.converted
        out_json_dict["in_models"] = self.in_models
        out_json_dict["annotation"] = self.annotation
        out_json_dict["reactions"] = self.reactions
        return out_json_dict

    def _to_json(self):  # GGE
        args_dict_out = self._args_to_dict()
        json_dict_out = {}
        keys_with_problems = ["reactions"]
        for key in args_dict_out.keys():
            if key not in keys_with_problems:
                json_dict_out[key] = args_dict_out[key]

        # Rewriting pointers to objects with just their names (ids?)
        json_dict_out["reactions"] = {}

        for c_key, c_val in args_dict_out["reactions"].items():
            if not c_val:  # dict or list is empty
                json_dict_out["reactions"][c_key] = c_val
                continue
            elif type(c_val) is list:
                json_dict_out["reactions"][c_key] = []
                for l in c_val:
                    json_dict_out["reactions"][c_key].append(f"{l._get_replace_tag()}")
                continue
            if c_key == "comparison":
                json_dict_out["reactions"][c_key] = {}
                continue
            elif type(c_val) is dict:
                raise Warning(f"Should not be a dict {c_val} for {c_key}")
                # json_dict_out["reactions"][c_key] = {}
                # for c_key2, c_val2 in c_val.items():
                #     json_dict_out["reactions"][c_key][c_key2] = []
                #     for el in c_val2:
                #         json_dict_out["reactions"][c_key][c_key2].append(
                #             f"{el._get_replace_tag()}"
                #         )
        return json.dumps(json_dict_out)

    def _get_replace_tag(self):  # GGE
        tag = ""
        if self.converted:
            tag = "~"
        else:
            tag = "%"  # unique tag for not converted, as they could have similar id as converted one
        return f"{tag}{self.id}"

    def _replace_tags_with_objects(self, all_objects_by_id: dict):  # GGE
        for c_key, c_val in self.reactions.items():
            if not c_val:  # dict or list is empty
                continue
            elif type(c_val) is list:
                for l in range(len(c_val)):
                    if c_val[l] in all_objects_by_id:
                        self.reactions[c_key][l] = all_objects_by_id[c_val[l]]
                    else:
                        self.reactions[c_key][l] = c_val[l]
                continue
            elif type(c_val) is dict:
                for c_key2, c_val2 in c_val.items():
                    for el in range(len(c_val2)):
                        if c_val2[el] in all_objects_by_id:
                            self.reactions[c_key][c_key2][el] = all_objects_by_id[
                                c_val2[el]
                            ]
                        else:
                            self.reactions[c_key][c_key2][el] = c_val2[el]
        return

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
    """Setting dictionaries for all genes selected for supermodel - self.converted and not selected - self.notconverted."""

    def __addNewGenes_conv(
        self,
        all_models_data: dict,
        gene_folder: PosixPath,
        do_old_genes: dict,
        do_mix_conv_notconv: bool,
    ):
        for model_id in list(all_models_data.keys()):
            blast_file = gene_folder / (model_id + "_blast.tsv")
            try:
                conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            except:
                warnings.warn(f"\nWarning! File {str(blast_file)} can't be opened.")
                if do_old_genes[model_id]:
                    print(f"\nOld gene will be used for {model_id}")
                    for gene in all_models_data[model_id]["preprocess_model"].genes:
                        if gene.id in self.assembly.keys():
                            self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                            getattr(self, model_id).update(
                                {gene.id: self.assembly.get(gene.id)}
                            )
                        else:
                            new_gene = NewGene(
                                False,
                                {
                                    "type": "NewGene",
                                    "args": {
                                        "new_id": gene.id,
                                        "old_id": gene.id,
                                        "source": model_id,
                                        "possible_sources": list(
                                            all_models_data.keys()
                                        ),
                                        "converted": False,
                                    },
                                },
                            )
                            self.assembly.update({gene.id: new_gene})
                            getattr(self, model_id).update({gene.id: new_gene})
                else:
                    print(f"No genes will be  used for {model_id}")
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
                                False,
                                {
                                    "type": "NewGene",
                                    "args": {
                                        "new_id": gene.id,
                                        "old_id": gene.id,
                                        "source": model_id,
                                        "possible_sources": list(
                                            all_models_data.keys()
                                        ),
                                        "converted": False,
                                    },
                                },
                            )
                            getattr(self, to_add).update({gene.id: new_gene})
                    elif type(attr.values[0]) != str:
                        if gene.id in getattr(self, to_add).keys():
                            getattr(self, to_add).get(gene.id)._updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                False,
                                {
                                    "type": "NewGene",
                                    "args": {
                                        "new_id": gene.id,
                                        "old_id": gene.id,
                                        "source": model_id,
                                        "possible_sources": list(
                                            all_models_data.keys()
                                        ),
                                        "converted": False,
                                    },
                                },
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
                                False,
                                {
                                    "type": "NewGene",
                                    "args": {
                                        "new_id": new_id,
                                        "old_id": gene.id,
                                        "source": model_id,
                                        "possible_sources": list(
                                            all_models_data.keys()
                                        ),
                                        "converted": True,
                                    },
                                },
                            )
                            self.assembly.update({new_id: new_gene})
                            getattr(self, model_id).update({new_id: new_gene})

    def __init__(
        self,
        is_loading: bool,
        args_dict,
        # structure: {"type": "SetofNewGenes", "args": {"all_models_data": , "gene_folder": ,  "do_mix_conv_notconv": }}
        #        all_models_data: dict,
        #        gene_folder,
        #        do_mix_conv_notconv: bool
    ):
        if args_dict["type"] != "SetofNewGenes":  # Houston, we have a problem
            raise ValueError(
                f"The args for the SetofNewGenes objects are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        if is_loading:
            print("loading SetofNewGenes...")
            args = json.loads(args)
            self.assembly = {}
            self.comparison = defaultdict(dict)
            self.notconverted = {}
            for a_key, a_val in args["assembly"].items():
                object_dict = json.loads(a_val)
                self.assembly.update(
                    {a_key: NewGene(True, {"type": "NewGene", "args": object_dict})}
                )

            # for a_key, a_val in args["comparison"].items():
            #     object_dict = json.loads(a_val)
            #     self.comparison.update(
            #         {a_key: NewGene(True, {"type": "NewGene", "args": object_dict})}
            #     )

            for a_key, a_val in args["notconverted"].items():
                object_dict = json.loads(a_val)
                self.notconverted.update(
                    {a_key: NewGene(True, {"type": "NewGene", "args": object_dict})}
                )
            return

        # GGE found it easier here to give pointers to each object (I'm becoming lasy)
        all_models_data = args["all_models_data"]
        gene_folder = args["gene_folder"]
        do_mix_conv_notconv = args["do_mix_conv_notconv"]
        do_old_genes = args["do_old_genes"]

        self.assembly = {}
        for source in list(all_models_data.keys()):
            setattr(self, source, {})
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        if gene_folder is not None:
            self.__addNewGenes_conv(
                all_models_data, gene_folder, do_old_genes, do_mix_conv_notconv
            )
        else:
            for model_id in list(all_models_data.keys()):
                if not do_old_genes[model_id]:
                    continue
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            False,
                            {
                                "type": "NewGene",
                                "args": {
                                    "new_id": gene.id,
                                    "old_id": gene.id,
                                    "source": model_id,
                                    "possible_sources": list(all_models_data.keys()),
                                    "converted": False,
                                },
                            },
                        )
                        self.assembly.update({gene.id: new_gene})
            for gene in self.assembly.values():
                for model_id in gene.in_models["models_list"]:
                    getattr(self, model_id).update({gene.id: gene})

    def _args_to_dict(self):  # GGE
        out_dict = {}
        out_dict["assembly"] = self.assembly
        out_dict["comparison"] = {}  # self.comparison
        out_dict["notconverted"] = self.notconverted
        return out_dict

    def _to_json(self):  # GGE
        args_dict_out = self._args_to_dict()
        json_dict_out = {"assembly": {}, "comparison": {}, "notconverted": {}}
        # rewrite each element as json, than no problem writing all args
        for a_key, a_val in args_dict_out["assembly"].items():
            json_dict_out["assembly"].update({a_key: a_val._to_json()})
        # for a_key, a_val in args_dict_out["comparison"].items():
        #     for a_key2, a_val2 in a_val.items():
        #         a_val.update({a_key2: a_val2._to_json()})
        for a_key, a_val in args_dict_out["notconverted"].items():
            json_dict_out["notconverted"].update({a_key: a_val._to_json()})
        return json.dumps(json_dict_out)

    def _generate_objects_dict(self):  # GGE
        all_objects_by_id = {}

        for a_key, a_val in self.assembly.items():
            if a_val._get_replace_tag() in all_objects_by_id:
                print("Already exists!!!! NewGen, assembly")
            all_objects_by_id[a_val._get_replace_tag()] = a_val

        #        comparison - can be ignored (no unique elements)

        for a_key, a_val in self.notconverted.items():
            if a_val._get_replace_tag() in all_objects_by_id:
                print(f"Already exists!!!! NewGen, notconverted: {a_val.id}")
            all_objects_by_id[a_val._get_replace_tag()] = a_val

        return all_objects_by_id

    def _replace_tilda_tags(self, all_objects_by_id: dict):  # GGE

        for a_key, a_val in self.assembly.items():
            a_val._replace_tags_with_objects(all_objects_by_id)

        for a_key, a_val in self.comparison.items():
            for a_key2, a_val2 in a_val.items():
                a_val2._replace_tags_with_objects(all_objects_by_id)

        for a_key, a_val in self.notconverted.items():
            a_val._replace_tags_with_objects(all_objects_by_id)

        return


class SuperModel:  # TODO REAL 30.08.23 add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values."""

    def __find_connections(
        self,
        connection_knowledge: KnowledgeConnectingOldNew,
        all_models_data: dict,
        do_mix: bool,
    ):
        for met in self.metabolites.assembly.values():
            met._find_reactions(connection_knowledge)
        for r in self.reactions.assembly.values():
            r._find_reactants_products(connection_knowledge, "reactants")
            r._find_reactants_products(connection_knowledge, "products")
            r._find_metabolites(connection_knowledge)
            g_to_add = r._find_gene_and_gpr(connection_knowledge)
            for model_id, gene_ids in g_to_add.items():
                for g_id in list(set(gene_ids)):
                    if g_id in self.genes.assembly.keys():
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
                g_to_add = r._find_gene_and_gpr(connection_knowledge, do_notconv=True)
                for model_id, gene_ids in g_to_add.items():
                    for g_id in list(set(gene_ids)):
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
                        met.charge.get(model_id).append(old_mets[0].charge)
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
        is_loading,
        args_dict,
        # structure: {"type": "SuperModel", "args":{"final_m_sel": , "final_m_not_sel": , "final_r_sel": , "final_r_not_sel": , "all_models_data": ,"additional_periplasmic_m": , "periplasmic_r": , "m_db_info": , "r_db_info": , "gene_folder": , "do_old_genes":, "do_mix_conv_notconv": , "and_as_solid": }}
        #        final_m_sel: dict, # Old input
        #        final_m_not_sel: dict,
        #        final_r_sel: dict,
        #        final_r_not_sel: dict,
        #        all_models_data: dict,
        #        additional_periplasmic_m: dict,
        #        periplasmic_r: dict,
        #        m_db_info: pd.core.frame.DataFrame,
        #        r_db_info: pd.core.frame.DataFrame,
        #        gene_folder,
        #        do_mix_conv_notconv: bool,
        #        and_as_solid: bool,
    ):
        if args_dict["type"] != "SuperModel":
            raise ValueError(
                f"The args for the SuperModel are from different object ({args_dict['type']})!"
            )

        args = args_dict["args"]

        if is_loading:
            print("loading supermodel...")
            self.sources = args["sources"]
            if "notes" in args.keys():
                self.notes = args["notes"]
            else:
                self.notes = {s: {} for s in self.sources}
            self.metabolites = SetofNewElements(
                True, {"type": "SetofNewElements", "args": args["metabolites"]}
            )
            self.reactions = SetofNewElements(
                True, {"type": "SetofNewElements", "args": args["reactions"]}
            )
            self.genes = SetofNewGenes(
                True, {"type": "SetofNewGenes", "args": args["genes"]}
            )

            # Gather all objects (by type) - id and pointer
            dict_objects_for_replace = {}
            dict_objects_for_replace.update(self.metabolites._generate_objects_dict())
            dict_objects_for_replace.update(self.reactions._generate_objects_dict())
            dict_objects_for_replace.update(self.genes._generate_objects_dict())

            # Replace all "replace_me_" tags (~) with object
            self.metabolites._replace_tilda_tags(dict_objects_for_replace)
            self.reactions._replace_tilda_tags(dict_objects_for_replace)
            self.genes._replace_tilda_tags(dict_objects_for_replace)

            # Add attributes for original models
            for source in self.sources:
                setattr(self.metabolites, source, {})
                setattr(self.reactions, source, {})
                setattr(self.genes, source, {})
            for new_id, new_obj in self.metabolites.assembly.items():
                for model_id in new_obj.in_models["models_list"]:
                    getattr(self.metabolites, model_id).update({new_id: new_obj})
            for new_id, new_obj in self.reactions.assembly.items():
                for model_id in new_obj.in_models["models_list"]:
                    getattr(self.reactions, model_id).update({new_id: new_obj})
            for new_id, new_obj in self.genes.assembly.items():
                for model_id in new_obj.in_models["models_list"]:
                    getattr(self.genes, model_id).update({new_id: new_obj})
            return

        # else: initial init

        # GGE: I am vary lasy now, dont want to change code below at all (so wrtiting all args as they were)
        final_m_sel = args["final_m_sel"]
        final_m_not_sel = args["final_m_not_sel"]
        final_r_sel = args["final_r_sel"]
        final_r_not_sel = args["final_r_not_sel"]
        all_models_data = args["all_models_data"]
        additional_periplasmic_m = args["additional_periplasmic_m"]
        periplasmic_r = args["periplasmic_r"]
        m_db_info = args["m_db_info"]
        r_db_info = args["r_db_info"]
        gene_folder = args["gene_folder"]
        do_old_genes = args["do_old_genes"]
        do_mix_conv_notconv = args["do_mix_conv_notconv"]
        and_as_solid = args["and_as_solid"]

        self.sources = list(all_models_data.keys())
        self.notes = {
            model_id: all_models_data[model_id]["preprocess_model"].notes
            for model_id in all_models_data.keys()
        }
        print("Creating metabolites for supermodel")
        self.metabolites = SetofNewElements(
            False,
            {
                "type": "SetofNewElements",
                "args": {
                    "element_type": "metabolites",
                    "selected": final_m_sel,
                    "not_selected": final_m_not_sel,
                    "model_ids": self.sources,
                    "db_info": m_db_info,
                    "do_mix_conv_notconv": do_mix_conv_notconv,
                },
            },
            additional_periplasmic_m,
        )
        print("Creating reactions for supermodel")
        self.reactions = SetofNewElements(
            False,
            {
                "type": "SetofNewElements",
                "args": {
                    "element_type": "reactions",
                    "selected": final_r_sel,
                    "not_selected": final_r_not_sel,
                    "model_ids": self.sources,
                    "db_info": r_db_info,
                    "do_mix_conv_notconv": do_mix_conv_notconv,
                },
            },
        )
        print("Creating genes for supermodel")
        self.genes = SetofNewGenes(
            False,
            {
                "type": "SetofNewGenes",
                "args": {
                    "all_models_data": all_models_data,
                    "gene_folder": gene_folder,
                    "do_old_genes": do_old_genes,
                    "do_mix_conv_notconv": do_mix_conv_notconv,
                },
            },
        )
        print("Connecting supermodel network")
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
                all_models_data,
                final_r_all,
                "reactions",
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
            do_old_genes,
        )
        self.__find_connections(
            connection_knowledge,
            all_models_data,
            do_mix_conv_notconv,
        )
        print("Finalizing supermodel attributes")
        self.__get_additional_attributes(
            self.sources, connection_knowledge, do_mix_conv_notconv
        )
        self.__runSwitchedMetabolites()
        self.__assemble_attributes(and_as_solid, do_mix_conv_notconv)

    def _args_to_dict(self):  # GGE
        out_json_dict = {}
        out_json_dict["sources"] = self.sources
        out_json_dict["notes"] = self.notes
        out_json_dict["metabolites"] = self.metabolites
        out_json_dict["reactions"] = self.reactions
        out_json_dict["genes"] = self.genes
        return out_json_dict

    def _to_json(self):  # GGE
        args_dict_out = self._args_to_dict()
        json_dict_out = {}
        json_dict_out["sources"] = args_dict_out["sources"]
        json_dict_out["notes"] = args_dict_out["notes"]
        json_dict_out["metabolites"] = args_dict_out["metabolites"]._to_json()
        json_dict_out["reactions"] = args_dict_out["reactions"]._to_json()
        json_dict_out["genes"] = args_dict_out["genes"]._to_json()
        return json.dumps(json_dict_out)

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
            coreN = getCore(self, number_of_model, operator.ge, and_as_solid)
            print(f"Results are saved in 'comparison' attribute as {coreN}")

    def exactly_in(self, number_of_model: int, and_as_solid=False):
        if (
            type(number_of_model) != int
            or number_of_model < 1
            or number_of_model > len(self.sources)
        ):
            raise ValueError("Number to check does not fit the number of models")
        else:
            coreN = getCore(self, number_of_model, operator.eq, and_as_solid)
            print(f"Results are saved in 'comparison' attribute as {coreN}")

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
                name = getDifference(self, yes, no, and_as_solid, short_name_len)
                print(f"Results are saved in 'comparison' attribute as {name}")

    def get_venn_segments(self, short_name_len=None, and_as_solid=False):
        """Getting metabolites and reactions networks for each Venn segment in Venn
        diagram."""
        if short_name_len is None:
            short_name_len = self.get_short_name_len()
        combinations = []
        for i in range(1, len(self.sources)):
            combinations.extend(itertools.combinations(self.sources, i))
        for combo in combinations:
            yes = sorted(list(combo))
            no = sorted((list(set(self.sources) - set(combo))))
            name = getDifference(self, yes, no, and_as_solid, short_name_len)
            print(f"Results are saved in 'comparison' attribute as {name}")

    def get_intersection(self, and_as_solid=False):
        coreN = getCore(self, len(self.sources), operator.ge, and_as_solid)
        print(f"Results are saved in 'comparison' attribute as {coreN}")

    def get_all_confidence_levels(self, and_as_solid=False):
        for i in range(len(self.sources), 1, -1):
            self.at_least_in(i, and_as_solid=and_as_solid)

    # def write_supermodel_to_pkl(self, output_name: str, recursion_limit=None):
    #     if not output_name.endswith(".pkl"):
    #         raise ValueError("Wrong extension of the file")
    #     if exists(output_name):
    #         raise ValueError("File already exist, change the name")
    #     else:
    #         # max_rec = 0x100000
    #         # resource.setrlimit(
    #         #     resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY]
    #         # )
    #         # sys.setrecursionlimit(max_rec)
    #         with open(output_name, "wb") as fh:
    #             dill.dump(self, fh)

    def write_supermodel_to_json(self, output_name: str):  # GGE
        if not output_name.endswith(".json"):
            raise ValueError("Wrong extension of the file")
        else:
            json_model_string = self._to_json()
            with open(output_name, "w") as outfile:
                outfile.write(json_model_string)


#
# def read_supermodel_from_pkl(input_name: str):
#     supermodel = dill.load(open(input_name, "rb"))
#     return supermodel


def read_supermodel_from_json(input_name: str):  # GGE
    with open(input_name, "r") as infile:
        load_dict = json.load(infile)
    supermodel = SuperModel(True, {"type": "SuperModel", "args": load_dict})
    return supermodel
