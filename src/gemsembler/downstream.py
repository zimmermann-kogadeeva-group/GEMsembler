import operator
import re
import warnings
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import dill
import matplotlib.patches as mpatches
import networkx as nx
import pandas as pd
from cobra.flux_analysis import pfba
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import leaves_list, linkage

plt.rcParams["svg.fonttype"] = "none"

import seaborn as sns

from .comparison import getCoreConnections, getCoreGPR
from .creation import SuperModel
from .drawing import (
    MET_NOT_INT_GLOBAL,
    custom_histplot,
    define_edge_features,
    define_node_features,
    draw_one_known_pathway,
    draw_one_synt_path,
    get_color_palette,
    get_pyvis_from_nx,
)

GLYCOLYSIS_GLOBAL = {
    "HEX1": [("glc__D_c", "g6p_c")],
    "PGI": [("g6p_c", "f6p_c")],
    "PFK": [("f6p_c", "fdp_c")],
    "FBA": [("fdp_c", "g3p_c")],
    "GAPD": [("g3p_c", "13dpg_c")],
    "PGK": [("13dpg_c", "3pg_c")],
    "PGM": [("3pg_c", "2pg_c")],
    "ENO": [("2pg_c", "pep_c")],
    "PYK": [("pep_c", "pyr_c")],
}
TCA_GLOBAL = {
    "PYK": [("pep_c", "pyr_c")],
    "PPC": [("pep_c", "oaa_c")],
    "PPCK": [("pep_c", "oaa_c")],
    "PEPCK_re": [("pep_c", "oaa_c")],
    "PC": [("pyr_c", "oaa_c")],
    "PDH": [("pyr_c", "accoa_c")],
    "PFL": [("pyr_c", "accoa_c")],
    "CS": [("accoa_c", "cit_c"), ("oaa_c", "cit_c")],
    "ACONT": [("cit_c", "icit_c")],
    "ACONTa": [("cit_c", "acon_C_c")],
    "ACONTb": [("acon_C_c", "icit_c")],
    "ICL": [("icit_c", "succ_c"), ("icit_c", "glx_c")],
    "ICDHyr": [("icit_c", "akg_c")],
    "ICDHx": [("icit_c", "akg_c")],
    "ICITRED": [("icit_c", "osuc_c")],
    "OSUCCL": [("osuc_c", "akg_c")],
    "AKGDH": [("akg_c", "succoa_c")],
    "OOR2r": [("akg_c", "succoa_c")],
    "AKGDa": [("akg_c", "sdhlam_c"), ("lpam_c", "sdhlam_c")],
    "AKGDb": [("sdhlam_c", "succoa_c"), ("sdhlam_c", "dhlam_c")],
    "PDHcr": [("dhlam_c", "lpam_c")],
    "SUCOAS": [("succoa_c", "succ_c")],
    "SUCDi": [("succ_c", "fum_c")],
    "FRD7": [("fum_c", "succ_c")],
    "FUM": [("fum_c", "mal__L_c")],
    "MALS": [("glx_c", "mal__L_c")],
    "MDH": [("mal__L_c", "oaa_c")],
    "MDH2": [("mal__L_c", "oaa_c")],
    "MDH3": [("mal__L_c", "oaa_c")],
}
PENTOSE_PHOSPHATE_PATHWAY_GLOBAL = {
    "G6PBDH": [("g6p_B_c", "6pgl_c")],
    "G6PDH2r": [("g6p_c", "6pgl_c")],
    "PGL": [("6pgl_c", "6pgc_c")],
    "GND": [("6pgc_c", "ru5p__D_c")],
    "RPE": [("ru5p__D_c", "xu5p__D_c")],
    "RPI": [("ru5p__D_c", "r5p_c")],
    "TKT1": [("xu5p__D_c", "g3p_c"), ("r5p_c", "s7p_c")],
    "TALA": [("g3p_c", "f6p_c"), ("s7p_c", "e4p_c")],
    "TKT2": [("e4p_c", "f6p_c"), ("xu5p__D_c", "g3p_c")],
    "PGI": [("f6p_c", "g6p_c")],
    "G6PI": [("g6p_c", "g6p_B_c")],
}
COFACTORS_GLOBAL = [
    "co2_c",
    "hco3_c",
    "pi_c",
    "ppi_c",
    "ACP_c",
    "atp_c",
    "adp_c",
    "amp_c",
    "nad_c",
    "nadh_c",
    "nadp_c",
    "nadph_c",
    "coa_c",
    "cmp_c",
    "cdp_c",
    "ctp_c",
    "gmp_c",
    "gdp_c",
    "gtp_c",
    "ump_c",
    "udp_c",
    "utp_c",
    "fadh2_c",
    "fad_c",
    "q8_c",
    "q8h2_c",
    "mqn8_c",
    "mql8_c",
    "mqn6_c",
    "mql6_c",
    "thf_c",
]


def write_metabolites_production_output(
    out_bp_production_tab,
    write_output_to_folder: str,
    plot_file_name=None,
    table_file_name=None,
    column_order=None,
    met_names=True,
    id_instead_long_name=20,
    cmap="cubehelix",
    metric="jaccard",
    method="single",
    linewidths=0.5,
    dpi=300,
    fig_height=None,
    fig_width=None,
    **kwargs,
):
    """Function to plot heatmap of produced or not produced metabolites.
    Parameters: table with production, output fold and
    attributes of seaborn clustermap function."""
    if table_file_name is None:
        out_bp_production_tab.to_csv(
            f"{write_output_to_folder}/all_metabolites_production.tsv",
            sep="\t",
            index=False,
        )
    else:
        out_bp_production_tab.to_csv(
            table_file_name,
            sep="\t",
            index=False,
        )
    if len(out_bp_production_tab.columns) <= 3:
        return []
    if met_names:
        if id_instead_long_name is not None:
            short_names = []
            for i in range(len(out_bp_production_tab["Metabolite names"])):
                if (
                    len(out_bp_production_tab["Metabolite names"].to_list()[i])
                    > id_instead_long_name
                ):
                    short_names.append(
                        out_bp_production_tab["Metabolites"].to_list()[i]
                    )
                else:
                    short_names.append(
                        out_bp_production_tab["Metabolite names"].to_list()[i]
                    )
            out_bp_production_tab.index = short_names
        else:
            out_bp_production_tab.index = out_bp_production_tab["Metabolite names"]
    else:
        out_bp_production_tab.index = out_bp_production_tab["Metabolites"]

    num_no_grow = out_bp_production_tab.drop(index="overall_growth", errors="ignore")
    num_no_grow = num_no_grow.astype(
        {
            mod_name: "float"
            for mod_name in num_no_grow.columns
            if mod_name not in ["Metabolites", "Metabolite names"]
        }
    )

    clust_df = num_no_grow.drop(columns=["Metabolites", "Metabolite names"])
    j = linkage(
        clust_df.values,
        method=method,
        metric=metric,
    )
    row_order = num_no_grow.index[leaves_list(j)]
    met_order = num_no_grow.loc[row_order]["Metabolites"].to_list()

    if column_order is None:
        jt = linkage(
            clust_df.values.transpose(),
            method=method,
            metric=metric,
        )
        column_order = clust_df.columns[leaves_list(jt)]

    df_plot = clust_df.loc[row_order, column_order]

    if fig_width is None:
        fig_width = 0.55 * len(column_order)
    if fig_height is None:
        fig_height = 0.17 * len(row_order)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cbar_ax = fig.add_axes([0.05, 0.05, 0.03, 0.1])
    fig.subplots_adjust(left=0.62, top=0.99, bottom=0.13, right=0.98)
    sns.heatmap(
        df_plot,
        cmap=sns.color_palette(cmap, 6),
        vmin=-0.5,
        vmax=5.5,
        cbar_ax=cbar_ax,
        cbar_kws=dict(orientation="vertical"),
        ax=ax,
        linewidth=linewidths,
        **kwargs,
    )
    cbar_ax.set_yticks(
        range(6),
        [
            "absent",
            "no synthesized",
            "not in biomass",
            "cofactor",
            "media",
            "synthesized",
        ],
        rotation=0,
        ha="left",
    )
    cbar_ax.set_title("Metabolite synthesis status")
    ax.tick_params(axis="both", which="major", pad=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    if plot_file_name is None:
        fig.savefig(f"{write_output_to_folder}/all_metabolites_production.png", dpi=dpi)
    else:
        fig.savefig(plot_file_name, dpi=dpi)
    return [fig, met_order]


def table_reactions_confidence(
    supermodel: SuperModel,
    output_name=None,
    pathway_r=None,
    path_r_dist=None,
    yes_range=1,
    no_range=1,
    genes=True,
    and_as_solid=False,
    add_original_models=True,
):
    """
    The function serves to organize information about reactions in a table format.
    Main columns of the table:
        Reaction - reaction ID
        R_confidence - reaction agreement score: in how many models the reaction is present
        Status - automatic status, based on agreement score; can be "yes" (rather present), "no" (rather not present), "q" (uncertain)
        R_models - in which models the reaction is present
        GPR_confidence - GPR agreement score: in how many models there is not empty common GPR part
        GPR_core - common GPR part from the most models
        GPR_assembly - united GPR from all models with the reaction present
    If pathway can be focused around one metabolite of interest, additional column:
        Distance from "met X" - how many reactions away is the reaction from metabolite of intererst
    If it is convenient to consider all input models (supormodel.sourcers) separately,
    additional columns for each input model ID
    """

    if pathway_r is None:
        pathway_r = list(supermodel.reactions.assembly.keys())
    output = {"Reaction": pathway_r, "R_confidence": [], "Status": [], "R_models": []}
    if path_r_dist is not None:
        if len(set(pathway_r) & set(list(path_r_dist.keys()))) != len(pathway_r):
            raise ValueError(
                f"Not all reactions ids from pathway are present"
                f" in the distance dictionary. "
                f"Please check pathway_r and path_r_dist input"
            )
        output.update({f"Distance from {path_r_dist['metabolite']}": []})
    if genes:
        output.update({"GPR_confidence": [], "GPR_core": [], "GPR_assembly": []})
    if add_original_models:
        for source in supermodel.sources:
            output.update({source: []})
    for r_id in pathway_r:
        if r_id in supermodel.reactions.assembly.keys():
            r = supermodel.reactions.assembly[r_id]
            output["R_confidence"].append(f"{r.in_models['models_amount']}")
            if r.in_models["models_amount"] >= len(supermodel.sources) - yes_range:
                output["Status"].append("yes")
            elif r.in_models["models_amount"] <= no_range:
                output["Status"].append("no")
            else:
                output["Status"].append("q")
            output["R_models"].append(" ".join(sorted(r.in_models["models_list"])))
            if genes:
                if not r.gene_reaction_rule["assembly"]:
                    output["GPR_confidence"].append("0")
                    output["GPR_core"].append("")
                    output["GPR_assembly"].append("")
                else:
                    for i in range(r.in_models["models_amount"], 0, -1):
                        gpr_core = getCoreGPR(
                            r.gene_reaction_rule,
                            i,
                            operator.ge,
                            r.in_models["models_list"],
                            and_as_solid,
                        )
                        if gpr_core:
                            output["GPR_confidence"].append(f"{i}")
                            output["GPR_core"].append(gpr_core[0])
                            output["GPR_assembly"].append(
                                r.gene_reaction_rule["assembly"][0]
                            )
                            break
            if add_original_models:
                for source in supermodel.sources:
                    if source in r.in_models["models_list"]:
                        output[source].append("+")
                    else:
                        output[source].append("-")
            if path_r_dist is not None:
                output[f"Distance from {path_r_dist['metabolite']}"].append(
                    path_r_dist[r_id]
                )
        else:
            output["R_confidence"].append("0")
            output["Status"].append("no")
            output["R_models"].append("")
            if genes:
                output["GPR_confidence"].append("0")
                output["GPR_core"].append("")
                output["GPR_assembly"].append("")
            if add_original_models:
                for source in supermodel.sources:
                    output[source].append("-")
            if path_r_dist is not None:
                output[f"Distance from {path_r_dist['metabolite']}"].append(
                    path_r_dist[r_id]
                )
    output_tb = pd.DataFrame(output)
    if output_name is not None:
        output_tb.to_csv(output_name, sep="\t", index=False)
    return output_tb


def pathway_of_interest(
    supermodel: SuperModel,
    pathway_r: dict or list,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
    add_original_models=True,
):
    if type(pathway_r) == list:
        if met_not_int is None:
            met_not_int = deepcopy(MET_NOT_INT_GLOBAL)
        pathway_r_dict = {}
        for r_id in pathway_r:
            pathway_r_dict[r_id] = []
            r = supermodel.reactions.assembly[r_id]
            for rea in r.reactants["assembly"]:
                for pro in r.products["assembly"]:
                    if (re.sub("_([cep])$", "", rea.id) not in met_not_int.keys()) and (
                        re.sub("_([cep])$", "", pro.id) not in met_not_int.keys()
                    ):
                        pathway_r_dict[r_id].append((rea.id, pro.id))
        pathway_r_list = pathway_r
    else:
        pathway_r_dict = pathway_r
        pathway_r_list = list(pathway_r.keys())
    if draw_pathway_to_file is not None:
        g = draw_one_known_pathway(
            supermodel,
            pathway_r_dict,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file is not None:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            pathway_r_list,
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
            add_original_models=add_original_models,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def biosynthesis_pathway_with_media_and_metabolite(
    supermodel: SuperModel,
    pathway: list,
    medium: list,
    met_to_synt: str,
    draw_pathway_to_file=None,
    write_pathway_table_to_file=None,
    calc_dist_from_synt_met=True,
    check_distance=5,
    yes_range=1,
    no_range=1,
    genes=True,
    and_as_solid=False,
    add_original_models=True,
    write_table=True,
    draw_met_not_int=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
):
    if draw_pathway_to_file is not None:
        g = draw_one_synt_path(
            supermodel,
            pathway,
            medium,
            [met_to_synt],
            draw_pathway_to_file,
            draw_met_not_int,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_pathway_table_to_file is not None:
        if calc_dist_from_synt_met:
            r_dist_dict = calc_dist_for_synt_path(
                pathway, met_to_synt, supermodel, check_distance
            )
            t = table_reactions_confidence(
                supermodel,
                write_pathway_table_to_file,
                pathway,
                r_dist_dict,
                yes_range,
                no_range,
                genes,
                and_as_solid,
                add_original_models,
            )
        else:
            t = table_reactions_confidence(
                supermodel,
                write_pathway_table_to_file,
                pathway,
                yes_range=yes_range,
                no_range=no_range,
                genes=genes,
                and_as_solid=and_as_solid,
                add_original_models=add_original_models,
            )
    if draw_pathway_to_file and write_pathway_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_pathway_table_to_file:
        return g
    elif not draw_pathway_to_file and write_pathway_table_to_file:
        return t


def glycolysis(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
    add_original_models=True,
):
    glycolysis = deepcopy(GLYCOLYSIS_GLOBAL)
    if draw_pathway_to_file is not None:
        g = draw_one_known_pathway(
            supermodel,
            glycolysis,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(glycolysis.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
            add_original_models=add_original_models,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def tca(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
    add_original_models=True,
):
    tca = deepcopy(TCA_GLOBAL)
    if draw_pathway_to_file:
        g = draw_one_known_pathway(
            supermodel,
            tca,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(tca.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
            add_original_models=add_original_models,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def pentose_phosphate(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
    add_original_models=True,
):
    pentose_phosphate_pathway = deepcopy(PENTOSE_PHOSPHATE_PATHWAY_GLOBAL)
    if draw_pathway_to_file:
        g = draw_one_known_pathway(
            supermodel,
            pentose_phosphate_pathway,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(pentose_phosphate_pathway.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
            add_original_models=add_original_models,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def biomass(
    supermodel: SuperModel,
    only_difference=False,
    draw_plot_to_file=None,
    write_table_to_file=None,
    directed=True,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    if (draw_plot_to_file is not None) and (not draw_plot_to_file.endswith(".html")):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    color_brewer = get_color_palette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    if write_table_to_file is not None:
        output = {
            "Metabolite": [],
            "Type": [],
            "Met_confidence": [],
            "Biomass_confidence": [],
            "Status": [],
        }
        for source in supermodel.sources:
            output.update({source: []})
    g = nx.DiGraph()
    biomass_r = supermodel.reactions.assembly.get("Biomass")
    colname_r = define_node_features(
        color_brewer,
        "single_path_r",
        biomass_r,
        n_letter,
    )
    g.add_node(
        colname_r[0],
        label=colname_r[1],
        shape="box",
        color=colname_r[2],
        title=colname_r[3],
    )
    core_rea = []
    if only_difference:
        core_rea = getCoreConnections(
            biomass_r.reactants,
            len(supermodel.sources),
            operator.ge,
            supermodel.sources,
        )
    for rea in biomass_r.reactants.get("assembly"):
        if rea not in core_rea:
            colname_rea = define_node_features(
                color_brewer, "metabolites", rea, n_letter
            )
            rea_edge = define_edge_features(
                color_brewer,
                "single_path_r",
                "metabolites",
                rea,
                colname_rea[1],
                biomass_r,
                colname_r[1],
                "reactants",
                n_letter,
            )
            g.add_node(
                colname_rea[0],
                shape="o",
                label=colname_rea[1],
                color=colname_rea[2],
                title=colname_rea[3],
            )
            for e in rea_edge:
                g.add_edge(
                    e[0],
                    e[1],
                    color=e[2],
                    font_color="black",
                    title=e[3],
                )
            if write_table_to_file is not None:
                output["Metabolite"].append(colname_rea[0])
                output["Type"].append("reactant")
                output["Met_confidence"].append(f"Core{rea.in_models['models_amount']}")
                cr = int(rea_edge[0][3].split(":\n")[0].split(" ")[-1])
                output["Biomass_confidence"].append(f"Core{cr}")
                if cr >= len(supermodel.sources) - yes_range:
                    output["Status"].append("yes")
                elif cr <= no_range:
                    output["Status"].append("no")
                else:
                    output["Status"].append("q")
                short_s_bp = rea_edge[0][3].split(":\n")[1]
                l = int(len(short_s_bp) / cr)
                sources_bp = [
                    short_s_bp[i : i + l] for i in range(0, len(short_s_bp), l)
                ]
                for source in supermodel.sources:
                    if source[:l] in sources_bp:
                        output[source].append("+")
                    else:
                        output[source].append("-")
    core_pro = []
    if only_difference:
        core_pro = getCoreConnections(
            biomass_r.products,
            len(supermodel.sources),
            operator.ge,
            supermodel.sources,
        )
    for pro in biomass_r.products.get("assembly"):
        if pro not in core_pro:
            colname_pro = define_node_features(
                color_brewer, "metabolites", pro, n_letter
            )
            pro_edge = define_edge_features(
                color_brewer,
                "single_path_r",
                "metabolites",
                pro,
                colname_pro[1],
                biomass_r,
                colname_r[1],
                "products",
                n_letter,
            )
            g.add_node(
                colname_pro[0],
                shape="o",
                label=colname_pro[1],
                color=colname_pro[2],
                title=colname_pro[3],
            )
            for e in pro_edge:
                g.add_edge(
                    e[0],
                    e[1],
                    color=e[2],
                    font_color="black",
                    title=e[3],
                )
            if write_table_to_file is not None:
                output["Metabolite"].append(colname_pro[0])
                output["Type"].append("product")
                output["Met_confidence"].append(f"Core{pro.in_models['models_amount']}")
                cp = int(pro_edge[0][3].split(":\n")[0].split(" ")[-1])
                output["Biomass_confidence"].append(f"Core{cp}")
                if cp >= len(supermodel.sources) - yes_range:
                    output["Status"].append("yes")
                elif cp <= no_range:
                    output["Status"].append("no")
                else:
                    output["Status"].append("q")
                short_s_bp = pro_edge[0][3].split(":\n")[1]
                l = int(len(short_s_bp) / cp)
                sources_bp = [
                    short_s_bp[i : i + l] for i in range(0, len(short_s_bp), l)
                ]
                for source in supermodel.sources:
                    if source[:l] in sources_bp:
                        output[source].append("+")
                    else:
                        output[source].append("-")
    if draw_plot_to_file:
        pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
        pyvis_graph.write_html(draw_plot_to_file, notebook=False)
    if write_table_to_file is not None:
        output_tb = pd.DataFrame(output)
        output_tb.to_csv(write_table_to_file, sep="\t", index=False)
        if draw_plot_to_file:
            return g, output_tb
        else:
            return output_tb
    elif draw_plot_to_file:
        return g


def get_met_neighborhood(
    supermodel: SuperModel,
    metabolite_id: str,
    neighborhood_dist=2,
    highly_connected_t=10,
    draw_neiborhood_to_file=None,
    draw_met_not_int=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
):
    if metabolite_id not in supermodel.metabolites.assembly.keys():
        raise ValueError(f"{metabolite_id} is not in the supermodel")
    if met_not_int is None:
        met_not_int = deepcopy(MET_NOT_INT_GLOBAL)
    all_r = set()
    all_g = set()
    med_high_connect = []
    all_m = {metabolite_id}
    dist_m = [metabolite_id]
    for i in range(neighborhood_dist):
        new_all_m = set()
        for m_id in dist_m:
            tmp_m = re.sub("_([cep])$", "", m_id)
            if (tmp_m in met_not_int.keys()) & (not draw_met_not_int):
                continue
            met = supermodel.metabolites.assembly.get(m_id)
            reactions = [r.id for r in met.reactions["assembly"]]
            if len(reactions) > highly_connected_t:
                reactions = []
                med_high_connect.append(re.sub("_([cep])$", "", met.id))
            all_r = all_r | set(reactions)
            for r_id in reactions:
                if r_id == "Biomass":
                    continue
                new_all_m = new_all_m | set(
                    [
                        m.id
                        for m in supermodel.reactions.assembly[r_id]
                        .metabolites["assembly"]
                        .keys()
                    ]
                )
                all_g = all_g | set(
                    [
                        g.id
                        for g in supermodel.reactions.assembly[r_id].genes["assembly"]
                    ]
                )
        dist_m = new_all_m
        all_m = all_m | new_all_m
    if draw_neiborhood_to_file is not None:
        graph = draw_one_synt_path(
            supermodel,
            list(all_r),
            med_high_connect,
            [metabolite_id],
            draw_neiborhood_to_file,
            draw_met_not_int,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
        return graph, all_r, all_g, all_m
    else:
        return all_r, all_g, all_m


def calc_dist_for_synt_path(
    pathway_rs: list,
    met_of_interest: str,
    supermodel: SuperModel,
    check_distance=5,
    highly_connected_t=50,
    draw_met_not_int=False,
):
    path_r_dist = {r: f">{check_distance}" for r in pathway_rs}
    path_r_dist.update({"metabolite": met_of_interest})
    for i in range(1, check_distance + 1):
        all_r, all_g, all_m = get_met_neighborhood(
            supermodel,
            met_of_interest,
            neighborhood_dist=i,
            highly_connected_t=highly_connected_t,
            draw_met_not_int=draw_met_not_int,
        )
        r_intersect = list(set(pathway_rs) & set(all_r))
        for ri in r_intersect:
            if path_r_dist[ri] == f">{check_distance}":
                path_r_dist[ri] = i
    return path_r_dist


def _fba_growth_met_production(
    models_to_analyse: dict,
    medium: dict,
    biomass_precursors=True,
    met_of_interest=None,
    biomass_r_id=None,
    flux_threshold=0.001,
):
    medium_no_comp = [re.sub("_([cep])$", "", m) for m in medium.keys()]
    if biomass_r_id is None:
        biomass_r_id = "Biomass"
    if met_of_interest is None and biomass_precursors:
        met_of_interest = []
        for m in models_to_analyse.values():
            met_of_interest = list(
                set(
                    met_of_interest
                    + [
                        react.id
                        for react in m.reactions.get_by_id(biomass_r_id).reactants
                    ]
                )
            )
    out_bp_production = {"Metabolites": ["overall_growth"] + met_of_interest}
    metname_of_interest = []
    for metid in met_of_interest:
        found = False
        for m in models_to_analyse.values():
            m_metid_all = [mm.id for mm in m.metabolites]
            if metid in m_metid_all:
                metname_of_interest.append(m.metabolites.get_by_id(metid).name)
                found = True
                break
        if not found:
            metname_of_interest.append(metid)
    out_bp_production.update(
        {"Metabolite names": ["overall_growth"] + metname_of_interest}
    )
    flux_res_out = defaultdict(dict)
    path_pfba_out = defaultdict(dict)
    stat_out = {}
    bp_synt = []
    for k, model in models_to_analyse.items():
        print(f"Running FBA block for {k} model")
        med_mod = {}
        for e in model.exchanges:
            if list(e.metabolites.keys())[0].id in medium.keys():
                med_mod.update({e.id: medium[list(e.metabolites.keys())[0].id]})
            else:
                med_mod.update({e.id: 0})
        old_medium = model.medium
        model.medium = med_mod
        res = model.optimize()
        flux_res_out[k] = {"overall_fba": res}
        model_data = [res.objective_value]
        model_stat = 0
        model_med_stat = 0
        all_met = [mm.id for mm in model.metabolites]
        all_r = [rr.id for rr in model.reactions]
        bp_model = [
            react.id for react in model.reactions.get_by_id(biomass_r_id).reactants
        ]
        for bp in met_of_interest:
            if bp not in all_met:
                model_data.append(0)
                continue
            demand_added = False
            if ("DM_" + bp) not in all_r:
                model.add_boundary(model.metabolites.get_by_id(bp), type="demand")
                demand_added = True
            model.objective = model.demands.get_by_id("DM_" + bp)
            res_bp = model.optimize()
            flux_res_out[k].update({bp + "_fba": res_bp})
            if res_bp.objective_value > flux_threshold:
                model_data.append(5)
                if re.sub("_([cep])$", "", bp) in medium_no_comp:
                    model_med_stat = model_med_stat + 1
                else:
                    model_stat = model_stat + 1
                if bp not in bp_synt:
                    bp_synt.append(bp)
                pfba_res = pfba(model)
                flux_res_out[k].update({bp + "_pfba": pfba_res})
                reactions = list(
                    pfba_res.to_frame()[
                        (pfba_res.to_frame()["fluxes"] > flux_threshold)
                        | (pfba_res.to_frame()["fluxes"] < -1 * flux_threshold)
                    ].index
                )
                path_pfba_out[k].update({bp: reactions})
            else:
                path_pfba_out[k].update({bp: ["No_path"]})
                if (bp not in bp_model) and biomass_precursors:
                    model_data.append(2)
                else:
                    model_data.append(1)
            if demand_added:
                model.remove_reactions(["DM_" + bp])
        model.medium = old_medium
        model.objective = model.reactions.get_by_id(biomass_r_id)
        out_bp_production.update({k: model_data})
        stat_out.update({k: model_stat})
        stat_out.update({"medium_" + k: model_med_stat})
    stat_out.update({"not_synthesized": len(set(met_of_interest) - set(bp_synt))})
    return out_bp_production, flux_res_out, path_pfba_out, stat_out


def _write_pfba_mq_results(
    path_pfba_mq_out: dict,
    supermodel: SuperModel,
    medium: list,
    output_folder: str,
    draw_pfba_mq_for_model=None,
    draw_met_not_int=False,
    draw_pfba_mq=True,
    table_pfba_mq=True,
    draw_confidence=True,
    output_folder_mq_paths_plots=None,
    output_folder_mq_paths_tables=None,
    confidence_plot=None,
    confidence_table=None,
    calc_r_dist=True,
    check_distance=5,
    met_order=None,
    met_rm_from_ptw=None,
    yes_range=1,
    no_range=1,
    genes=True,
    and_as_solid=False,
    add_original_models=True,
    dpi=300,
    fig_height=None,
    **kwargs,
):
    if (draw_pfba_mq is True) and (output_folder_mq_paths_plots is None):
        if output_folder_mq_paths_tables is not None:
            output_folder_mq_paths_plots = output_folder_mq_paths_tables
        else:
            output_folder_mq_paths_plots = output_folder
    if (table_pfba_mq is True) and (output_folder_mq_paths_tables is None):
        if output_folder_mq_paths_plots is not None:
            output_folder_mq_paths_tables = output_folder_mq_paths_plots
        else:
            output_folder_mq_paths_tables = output_folder
    if draw_pfba_mq_for_model is None:
        met_model = {}
        all_met = []
        for mod, res in path_pfba_mq_out.items():
            for met in res.keys():
                if met not in all_met:
                    all_met.append(met)
                if res[met] != ["No_path"]:
                    if met in met_model.keys():
                        if mod.startswith("core") and (
                            (not met_model[met][0].startswith("core"))
                            or (
                                int(mod.removeprefix("core"))
                                > int(met_model[met][0].removeprefix("core"))
                            )
                        ):
                            met_model[met] = [mod]
                        elif mod == "assembly" and (
                            not met_model[met][0].startswith("core")
                        ):
                            met_model[met] = [mod]
                        elif (not met_model[met][0].startswith("core")) and (
                            met_model[met][0] != "assembly"
                        ):
                            met_model[met].append(mod)
                    else:
                        met_model[met] = [mod]
        pd.DataFrame(
            met_model.items(), columns=["Metabolite", "Most confident path"]
        ).to_csv(
            f"{output_folder}/bp_most_confident_met_production.tsv",
            index=False,
            sep="\t",
        )
    else:
        met_model = {}
        all_met = []
        for met in path_pfba_mq_out[draw_pfba_mq_for_model].keys():
            if met not in all_met:
                all_met.append(met)
            if path_pfba_mq_out[draw_pfba_mq_for_model][met] != ["No_path"]:
                met_model[met] = [draw_pfba_mq_for_model]
    if draw_confidence:
        confidence_paths = {
            "ID": [],
            "Reactions/GPRs": [],
            "Metabolite synthesis": [],
            "Confidence": [],
        }
        if met_order is not None:
            all_met = met_order
    if met_rm_from_ptw is not None:
        all_met = [m for m in all_met if m not in met_rm_from_ptw]
    for m in all_met:
        if m not in met_model.keys():
            if draw_confidence:
                confidence_paths["ID"].append("No_r")
                confidence_paths["Reactions/GPRs"].append("Reactions")
                confidence_paths["Metabolite synthesis"].append(m)
                confidence_paths["Confidence"].append(0)
                confidence_paths["ID"].append("No_gpr")
                confidence_paths["Reactions/GPRs"].append("GPRs")
                confidence_paths["Metabolite synthesis"].append(m)
                confidence_paths["Confidence"].append(0)
            continue
        for model_id in met_model[m]:
            v = [vv for vv in path_pfba_mq_out[model_id][m] if not vv.startswith("DM_")]
            if draw_pfba_mq:
                g = draw_one_synt_path(
                    supermodel,
                    v,
                    medium,
                    [m],
                    f"{output_folder_mq_paths_plots}/{m}_{model_id}.html",
                    draw_met_not_int=draw_met_not_int,
                )
            if table_pfba_mq:
                if calc_r_dist:
                    r_dist_dict = calc_dist_for_synt_path(
                        v, m, supermodel, check_distance
                    )
                else:
                    r_dist_dict = None
                t = table_reactions_confidence(
                    supermodel,
                    f"{output_folder_mq_paths_tables}/{m}_{model_id}.tsv",
                    v,
                    r_dist_dict,
                    yes_range,
                    no_range,
                    genes,
                    and_as_solid,
                    add_original_models,
                )

            if draw_confidence:
                for vr in v:
                    confidence_paths["ID"].append(vr)
                    confidence_paths["Reactions/GPRs"].append("Reactions")
                    confidence_paths["Metabolite synthesis"].append(m)
                    confidence_paths["Confidence"].append(
                        supermodel.reactions.assembly[vr].in_models["models_amount"]
                    )
                    if not supermodel.reactions.assembly[vr].gene_reaction_rule[
                        "assembly"
                    ]:
                        confidence_paths["ID"].append("No_gpr")
                        confidence_paths["Reactions/GPRs"].append("GPRs")
                        confidence_paths["Metabolite synthesis"].append(m)
                        confidence_paths["Confidence"].append(0)
                    else:
                        for i in range(
                            supermodel.reactions.assembly[vr].in_models[
                                "models_amount"
                            ],
                            0,
                            -1,
                        ):
                            gpr_core = getCoreGPR(
                                supermodel.reactions.assembly[vr].gene_reaction_rule,
                                i,
                                operator.ge,
                                supermodel.reactions.assembly[vr].in_models[
                                    "models_list"
                                ],
                                False,
                            )
                            if gpr_core:
                                confidence_paths["ID"].append(gpr_core[0])
                                confidence_paths["Reactions/GPRs"].append("GPRs")
                                confidence_paths["Metabolite synthesis"].append(
                                    m.removesuffix("_path_pfba")
                                )
                                confidence_paths["Confidence"].append(i)
                                break
    if draw_confidence:
        confidence_paths_tab = pd.DataFrame(confidence_paths)
        if confidence_table is not None:
            confidence_paths_tab.to_csv(
                f"{confidence_table}",
                index=False,
                sep="\t",
            )
        else:
            confidence_paths_tab.to_csv(
                f"{output_folder}/confidence_met_production_paths.tsv",
                index=False,
                sep="\t",
            )
        if fig_height is None:
            fig_height = 0.17 * len(all_met)
        g = (
            sns.FacetGrid(
                data=confidence_paths_tab,
                col="Reactions/GPRs",
                height=fig_height,
                aspect=0.5,
                legend_out=True,
            )
            .map_dataframe(
                custom_histplot,
                y="Metabolite synthesis",
                hue="Confidence",
                **kwargs,
            )
            .set_titles(col_template="{col_name}")
            .set(
                ylim=(
                    confidence_paths_tab.get("Metabolite synthesis").unique().shape[0]
                    - 0.5,
                    -0.5,
                )
            )
        )
        uniq_r = list(
            confidence_paths_tab[confidence_paths_tab["Reactions/GPRs"] == "Reactions"][
                "Confidence"
            ]
            .sort_values(ascending=False)
            .unique()
        )
        uniq_gpr = list(
            confidence_paths_tab[confidence_paths_tab["Reactions/GPRs"] == "GPRs"][
                "Confidence"
            ]
            .sort_values(ascending=False)
            .unique()
        )
        g.axes[0, 1].legend(
            handles=[x[0] for x in g.axes[0, 0].containers]
            + [mpatches.Patch(color="none", label="")]
            * abs(len(uniq_r) - len(uniq_gpr))
            + [x[0] for x in g.axes[0, 1].containers],
            labels=uniq_r + [""] * abs(len(uniq_r) - len(uniq_gpr)) + uniq_gpr,
            ncol=2,
            loc="center right",
            bbox_to_anchor=(1.4, 0.1),
            title="Confidence level\nReactions   GPRs",
        )
        g.figure.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.97)
        if confidence_plot is not None:
            g.savefig(confidence_plot, dpi=dpi)
        else:
            g.savefig(f"{output_folder}/confidence_met_production_paths.png", dpi=dpi)
        return [g]
    else:
        return []


def run_growth_full_flux_analysis(
    models_to_analyse: dict,
    medium: dict,
    supermodel: SuperModel,
    output_folder: str,
    biomass_precursors=True,
    metabolites_of_interest=None,
    production_plot=None,
    production_table=None,
    stat_file=None,
    model_to_further_analyse=None,
    draw_pfba=True,
    table_pfba=True,
    draw_confidence=True,
    output_folder_mq_paths_plots=None,
    output_folder_mq_paths_tables=None,
    confidence_plot=None,
    confidence_table=None,
    calc_dist_path=True,
    check_distance=5,
    draw_met_not_int=False,
    biomass_r_id=None,
    met_names=True,
    column_order=None,
    met_rm_from_ptw=None,
    id_instead_long_name=20,
    dpi=300,
    fig_height=None,
    fig_width=None,
    flux_threshold=0.001,
    **kwargs,
):
    """"""
    (
        out_bp_production,
        flux_res_out,
        path_pfba_out,
        stat_out,
    ) = _fba_growth_met_production(
        models_to_analyse,
        medium,
        biomass_precursors,
        metabolites_of_interest,
        biomass_r_id,
        flux_threshold=flux_threshold,
    )
    out_bp_production_tab = pd.DataFrame(out_bp_production)
    stat_out_tab = pd.DataFrame(
        stat_out.items(),
        columns=["Metabolites confidence production", "Metabolites amount"],
    )
    production_plots = write_metabolites_production_output(
        out_bp_production_tab,
        output_folder,
        production_plot,
        production_table,
        met_names=met_names,
        column_order=column_order,
        id_instead_long_name=id_instead_long_name,
        dpi=dpi,
        fig_height=fig_height,
        fig_width=fig_width,
        **kwargs,
    )
    if not production_plots:
        met_order = out_bp_production_tab.drop(index="overall_growth", errors="ignore")[
            "Metabolites"
        ].to_list()
    else:
        met_order = production_plots[1]
    if stat_file is None:
        stat_file = f"{output_folder}/production_confidence_stat.tsv"
    stat_out_tab.to_csv(stat_file, sep="\t", index=False)
    if draw_pfba or table_pfba or draw_confidence:
        if model_to_further_analyse is None:
            warnings.warn(
                "Models for which to draw pfba of biomass precursors is not provided. "
                "So for each biomass precursor the model with "
                "the highest confidence level will be used."
            )
        production_plots = production_plots + _write_pfba_mq_results(
            path_pfba_out,
            supermodel,
            list(medium.keys()),
            output_folder,
            model_to_further_analyse,
            draw_met_not_int,
            draw_pfba,
            table_pfba,
            draw_confidence,
            output_folder_mq_paths_plots,
            output_folder_mq_paths_tables,
            confidence_plot,
            confidence_table,
            calc_dist_path,
            check_distance,
            met_order,
            met_rm_from_ptw,
            dpi,
            fig_height,
            **kwargs,
        )
    return (
        out_bp_production_tab,
        production_plots,
        flux_res_out,
        path_pfba_out,
        stat_out_tab,
    )


def run_metquest_results_analysis(
    folder_with_mq_res_folders: str,
    model_list: list,
    metabolites_ids: list,
    medium: list,
    supermodel: SuperModel,
    output_folder: str,
    cofactors=None,
    production_plot=None,
    production_table=None,
    stat_file=None,
    model_to_further_analyse=None,
    draw_mq_path=False,
    table_mq_path=False,
    draw_confidence=False,
    output_folder_mq_paths_plots=None,
    output_folder_mq_paths_tables=None,
    confidence_plot=None,
    confidence_table=None,
    calc_dist_path=True,
    check_distance=5,
    draw_met_not_int=False,
    check_in_biomass_precursors=False,
    met_names=True,
    column_order=None,
    met_rm_from_ptw=None,
    id_instead_long_name=20,
    dpi=300,
    fig_height=None,
    fig_width=None,
    **kwargs,
):
    if output_folder_mq_paths_plots is not None:
        draw_mq_path = True
    if output_folder_mq_paths_tables is not None:
        table_mq_path = True
    if cofactors is None:
        cofactors = deepcopy(COFACTORS_GLOBAL)
    if check_in_biomass_precursors:
        all_models_bp = {}
        for model in model_list:
            if model in supermodel.sources + ["assembly"]:
                all_models_bp[model] = [
                    m.id
                    for m in supermodel.reactions.assembly["Biomass"].reactants.get(
                        model, []
                    )
                ]
            else:
                all_models_bp[model] = [
                    m.id
                    for m in supermodel.reactions.assembly["Biomass"]
                    .reactants["comparison"]
                    .get(model, [])
                ]
    metquest_all_res_paths = {}
    stat_out = {}
    met_synt = []
    met_med = []
    met_cof = []
    folder_with_out_folders = Path(folder_with_mq_res_folders)
    for mod in model_list:
        for f in folder_with_out_folders.iterdir():
            if f.is_dir() and f.name.startswith(mod):
                metquest_all_res_paths[mod] = f / "shortest_paths.pkl"
    metabolites_names = []
    for mid in metabolites_ids:
        if mid in supermodel.metabolites.assembly.keys():
            metabolites_names.append(supermodel.metabolites.assembly[mid].name)
        else:
            metabolites_names.append(mid)
    synthes_out = {
        "Metabolites": metabolites_ids,
        "Metabolite names": metabolites_names,
    } | {mod: [] for mod in model_list}
    met_interest_mq_paths = {}
    for mod in model_list:
        met_interest_mq_paths[mod] = {}
        model_stat = 0
        mq_data = dill.load(open(metquest_all_res_paths[mod], "rb"))
        for bp in metabolites_ids:
            if bp not in mq_data["paths"].keys():
                to_append = 0
            else:
                if bp in medium:
                    to_append = 4
                    if bp not in met_med:
                        met_med.append(bp)
                elif bp in cofactors:
                    to_append = 3
                    if bp not in met_cof:
                        met_cof.append(bp)
                elif mq_data["paths"][bp]:
                    to_append = 5
                    model_stat = model_stat + 1
                    if bp not in met_synt:
                        met_synt.append(bp)
                elif check_in_biomass_precursors:
                    if bp not in all_models_bp[mod]:
                        to_append = 2
                    else:
                        to_append = 1
                else:
                    to_append = 1
            synthes_out[mod].append(to_append)
            if to_append == 5:
                allalt = set()
                for altp_len in mq_data["paths"][bp].values():
                    for altp in altp_len:
                        allalt = allalt | altp
                met_interest_mq_paths[mod].update({bp: list(allalt)})
            else:
                met_interest_mq_paths[mod].update({bp: ["No_path"]})
        stat_out.update({mod: model_stat})
    synthes_tab_out = pd.DataFrame(synthes_out)
    stat_out.update(
        {
            "Metabolites in medium": len(met_med),
            "Metabolites in cofactors": len(met_cof),
            "Metabolites not synthesized": len(
                set(metabolites_ids) - set(met_synt) - set(met_med) - set(met_cof)
            ),
        }
    )
    stat_out_tab = pd.DataFrame(
        stat_out.items(),
        columns=["Metabolites confidence production", "Metabolites amount"],
    )
    production_plots = write_metabolites_production_output(
        synthes_tab_out,
        output_folder,
        production_plot,
        production_table,
        met_names=met_names,
        column_order=column_order,
        id_instead_long_name=id_instead_long_name,
        dpi=dpi,
        fig_height=fig_height,
        fig_width=fig_width,
        **kwargs,
    )
    if not production_plots:
        met_order = synthes_tab_out["Metabolites"].to_list()
        if not column_order:
            column_order = list(synthes_tab_out.columns)
    else:
        met_order = production_plots[1]
        column_order = production_plots[2]
    if stat_file is None:
        stat_file = f"{output_folder}/production_confidence_stat.tsv"
    stat_out_tab.to_csv(stat_file, sep="\t", index=False)
    if draw_mq_path or table_mq_path or draw_confidence:
        production_plots = production_plots + _write_pfba_mq_results(
            met_interest_mq_paths,
            supermodel,
            medium,
            output_folder,
            model_to_further_analyse,
            draw_met_not_int,
            draw_mq_path,
            table_mq_path,
            draw_confidence,
            output_folder_mq_paths_plots,
            output_folder_mq_paths_tables,
            confidence_plot,
            confidence_table,
            calc_dist_path,
            check_distance,
            met_order,
            met_rm_from_ptw,
            dpi,
            fig_height,
            **kwargs,
        )
    return synthes_tab_out, production_plots, stat_out_tab, met_interest_mq_paths
