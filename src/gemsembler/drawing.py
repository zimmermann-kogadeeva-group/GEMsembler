import math
import operator
import re
import networkx as nx
import pandas as pd
import seaborn as sns
from .creation import NewElement, SuperModel
from .comparison import getCoreGPR, getCoreConnections
from pyvis.network import Network


def get_color_palette(sources_number: int) -> dict:
    up = math.ceil(sources_number / 10)
    down = math.ceil(sources_number / 10)
    colors = sns.color_palette("bright").as_hex()
    out_colors = {
        "metabolites": sns.light_palette(
            colors[4], sources_number + up + down
        ).as_hex()[up:],
        "intersection_r": sns.light_palette(
            colors[3], sources_number + up + down
        ).as_hex()[up:],
        "single_path_r": sns.light_palette(
            colors[1], sources_number + up + down
        ).as_hex()[up:-down],
        "genes": sns.light_palette(colors[8], sources_number + up + down).as_hex()[
            up:-down
        ],
        "metConnect": sns.light_palette(colors[0], sources_number + up + down).as_hex()[
            up:-down
        ],
        "reactionsConnect": sns.light_palette(
            colors[9], sources_number + up + down
        ).as_hex()[up:-down],
        "interest": sns.light_palette(colors[2], sources_number + up + down).as_hex()[
            up:-down
        ],
        "notFound": sns.light_palette(colors[7], sources_number + up + down).as_hex()[
            up:-down
        ],
    }
    return out_colors


def define_node_features(
    colordata: dict,
    pallitra: str,
    object: NewElement,
    n_letter: int,
    gene=False,
    and_as_solid=False,
    node_id=None,
    label=None,
) -> [str]:
    if not gene:
        col = colordata.get(pallitra)[object.in_models["models_amount"] - 1]
        short_name = "".join(
            [tmp[:n_letter] for tmp in object.in_models["models_list"]]
        )
        if node_id is None:
            node_id = object.id
        if label is None:
            label = object.id
        title = f"Confirmed by {object.in_models['models_amount']}:\n{short_name}"
    else:
        if not object.gene_reaction_rule["assembly"]:
            node_id = f"{object.id} NO GPR"
            label = "NO GPR"
            col = colordata.get("notFound")[object.in_models["models_amount"] - 1]
            title = "No genes found"
        else:
            for i in range(object.in_models["models_amount"], 0, -1):
                gpr = getCoreGPR(
                    object.gene_reaction_rule,
                    i,
                    operator.ge,
                    object.in_models["models_list"],
                    and_as_solid,
                )
                if gpr:
                    node_id = f"{object.id} GPR core {i}"
                    label = f"GPR core {i}"
                    col = colordata.get(pallitra)[i - 1]
                    title = (
                        f"GPR core {i}:\n{gpr[0]}\n"
                        f"GPR assembly:\n{object.gene_reaction_rule['assembly'][0]}"
                    )
                    break
    return [node_id, label, col, title]


def define_edge_features(
    colordata: dict,
    pallitra_forward: str,
    pallitra_reversed: str,
    mr: NewElement,
    colname_m: str,
    r: NewElement,
    colname_r: str,
    connections: str,
    n_letter: int,
) -> [[str, str, str, str]]:
    n = 0
    short_name = ""
    if r.lower_bound.get("assembly")[0] >= 0:
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp):
                n = n + 1
                short_name = short_name + tmp[:n_letter]
        if connections == "reactants":
            col = colordata.get(pallitra_forward)[n - 1]
            output = [[colname_m, colname_r, col, f"Confirmed by {n}:\n{short_name}"]]
        elif connections == "products":
            col = colordata.get(pallitra_reversed)[n - 1]
            output = [[colname_r, colname_m, col, f"Confirmed by {n}:\n{short_name}"]]
    elif r.upper_bound.get("assembly")[0] <= 0:
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp):
                n = n + 1
                short_name = short_name + tmp[:n_letter]
        if connections == "reactants":
            col = colordata.get(pallitra_reversed)[n - 1]
            output = [[colname_r, colname_m, col, f"Confirmed by {n}:\n{short_name}"]]
        elif connections == "products":
            col = colordata.get(pallitra_forward)[n - 1]
            output = [[colname_m, colname_r, col, f"Confirmed by {n}:\n{short_name}"]]
    else:
        f = 0
        b = 0
        short_name_f = ""
        short_name_b = ""
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp) and r.upper_bound.get(tmp)[0] > 0:
                f = f + 1
                short_name_f = short_name_f + tmp[:n_letter]
            if mr in getattr(r, connections).get(tmp) and r.lower_bound.get(tmp)[0] < 0:
                b = b + 1
                short_name_b = short_name_b + tmp[:n_letter]
        if connections == "reactants":
            col_f = colordata.get(pallitra_forward)[f - 1]
            col_b = colordata.get(pallitra_forward)[b - 1]
            output = [
                [colname_m, colname_r, col_f, f"Confirmed by {f}:\n{short_name_f}"],
                [colname_r, colname_m, col_b, f"Confirmed by {b}:\n{short_name_b}"],
            ]
        elif connections == "products":
            col_f = colordata.get(pallitra_reversed)[f - 1]
            col_b = colordata.get(pallitra_reversed)[b - 1]
            output = [
                [colname_r, colname_m, col_f, f"Confirmed by {f}:\n{short_name_f}"],
                [colname_m, colname_r, col_b, f"Confirmed by {b}:\n{short_name_b}"],
            ]
    return output


def get_pyvis_from_nx(g, directed, size, width=1920, height=1080):
    pyvis_graph = Network(
        width=f"{width}px",
        height=f"{height}px",
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    for n in pyvis_graph.nodes:
        n["size"] = size
        n["font"] = {"size": size}
    return pyvis_graph


def draw_one_known_pathway(
    supermodel: SuperModel,
    pathway: dict,
    output_name: str,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
    size=25,
):
    if not output_name.endswith(".html"):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    if met_not_int is None:
        met_not_int = {
            "h": 0,
            "h2o": 0,
            "h2": 0,
            "oh1": 0,
            "o2": 0,
            "co2": 0,
            "coa": 0,
            "ppi": 0,
            "pi": 0,
            "amp": 0,
            "adp": 0,
            "atp": 0,
            "cmp": 0,
            "cdp": 0,
            "ctp": 0,
            "gmp": 0,
            "gdp": 0,
            "gtp": 0,
            "ump": 0,
            "udp": 0,
            "utp": 0,
            "nad": 0,
            "nadh": 0,
            "nadp": 0,
            "nadph": 0,
            "dadp": 0,
            "damp": 0,
            "nh3": 0,
            "nh4": 0,
            "fadh2": 0,
            "fad": 0,
            "ac": 0,
            "accoa": 0,
            "h2s": 0,
            "HC00250": 0,
        }
    color_brewer = get_color_palette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    g = nx.DiGraph()
    path_met = []
    for value in pathway.values():
        for v in value:
            path_met.append(v[0])
            path_met.append(v[1])
    path_met = list(set(path_met))
    for r_id in pathway.keys():
        if r_id in supermodel.reactions.assembly.keys():
            r = supermodel.reactions.assembly.get(r_id)
            colname_r = define_node_features(
                color_brewer,
                "single_path_r",
                r,
                n_letter,
            )
            g.add_node(
                colname_r[0],
                label=colname_r[1],
                shape="box",
                color=colname_r[2],
                title=colname_r[3],
            )
            for rea in r.reactants.get("assembly"):
                tmp_rea = re.sub("_([cep])$", "", rea.id)
                if rea.id in path_met:
                    colname_rea = define_node_features(
                        color_brewer,
                        "metabolites",
                        rea,
                        n_letter,
                    )
                    rea_edge = define_edge_features(
                        color_brewer,
                        "single_path_r",
                        "metabolites",
                        rea,
                        colname_rea[1],
                        r,
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
                        )
                elif additional_met:
                    if tmp_rea not in met_not_int.keys():
                        colname_rea = define_node_features(
                            color_brewer,
                            "metConnect",
                            rea,
                            n_letter,
                        )
                        rea_edge = define_edge_features(
                            color_brewer,
                            "single_path_r",
                            "metabolites",
                            rea,
                            colname_rea[1],
                            r,
                            colname_r[1],
                            "reactants",
                            n_letter,
                        )
                    else:
                        colname_rea = define_node_features(
                            color_brewer,
                            "metConnect",
                            rea,
                            n_letter,
                            node_id=f"{tmp_rea} {met_not_int[tmp_rea]}",
                            label=tmp_rea,
                        )
                        rea_edge = define_edge_features(
                            color_brewer,
                            "single_path_r",
                            "metabolites",
                            rea,
                            f"{tmp_rea} {met_not_int[tmp_rea]}",
                            r,
                            colname_r[1],
                            "reactants",
                            n_letter,
                        )
                        met_not_int[tmp_rea] = met_not_int[tmp_rea] + 1
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
                        )
            for pro in r.products.get("assembly"):
                tmp_pro = re.sub("_([cep])$", "", pro.id)
                if pro.id in path_met:
                    colname_pro = define_node_features(
                        color_brewer,
                        "metabolites",
                        pro,
                        n_letter,
                    )
                    pro_edge = define_edge_features(
                        color_brewer,
                        "single_path_r",
                        "metabolites",
                        pro,
                        colname_pro[1],
                        r,
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
                        )
                elif additional_met:
                    if tmp_pro not in met_not_int.keys():
                        colname_pro = define_node_features(
                            color_brewer,
                            "metConnect",
                            pro,
                            n_letter,
                        )
                        pro_edge = define_edge_features(
                            color_brewer,
                            "single_path_r",
                            "metabolites",
                            pro,
                            colname_pro[1],
                            r,
                            colname_r[1],
                            "products",
                            n_letter,
                        )
                    else:
                        colname_pro = define_node_features(
                            color_brewer,
                            "metConnect",
                            pro,
                            n_letter,
                            node_id=f"{tmp_pro} {met_not_int[tmp_pro]}",
                            label=tmp_pro,
                        )
                        pro_edge = define_edge_features(
                            color_brewer,
                            "single_path_r",
                            "metabolites",
                            pro,
                            f"{tmp_pro} {met_not_int[tmp_pro]}",
                            r,
                            colname_r[1],
                            "products",
                            n_letter,
                        )
                        met_not_int[tmp_pro] = met_not_int[tmp_pro] + 1
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
                        )
            if genes:
                g_colname = define_node_features(
                    color_brewer,
                    "genes",
                    r,
                    n_letter,
                    gene=True,
                    and_as_solid=and_as_solid,
                )
                g.add_node(
                    g_colname[0],
                    shape="box",
                    label=g_colname[1],
                    color=g_colname[2],
                    title=g_colname[3],
                )
                g.add_edge(
                    r_id,
                    g_colname[0],
                    color=g_colname[2],
                    font_color="black",
                )
        else:
            not_f_col = color_brewer["notFound"][
                math.floor(len(supermodel.sources) / 2)
            ]
            g.add_node(
                r_id,
                label=r_id,
                shape="box",
                color=not_f_col,
                title="Confirmed by:\nNone",
            )
            for pair in pathway.get(r_id):
                if pair[0] in supermodel.metabolites.assembly.keys():
                    colname_rea = define_node_features(
                        color_brewer,
                        "metabolites",
                        supermodel.metabolites.assembly[pair[0]],
                        n_letter,
                    )
                else:
                    colname_rea = [pair[0], pair[0], not_f_col, "Confirmed by:\nNone"]
                if pair[1] in supermodel.metabolites.assembly.keys():
                    colname_pro = define_node_features(
                        color_brewer,
                        "metabolites",
                        supermodel.metabolites.assembly[pair[1]],
                        n_letter,
                    )
                else:
                    colname_pro = [pair[1], pair[1], not_f_col, "Confirmed by:\nNone"]
                g.add_node(
                    colname_rea[0],
                    label=colname_rea[1],
                    shape="o",
                    color=colname_rea[2],
                    title=colname_rea[3],
                )
                g.add_node(
                    colname_pro[0],
                    label=colname_pro[1],
                    shape="o",
                    color=colname_pro[2],
                    title=colname_pro[3],
                )
                g.add_edge(
                    colname_rea[1],
                    r_id,
                    color=not_f_col,
                    font_color="black",
                )
                g.add_edge(
                    r_id,
                    colname_pro[1],
                    color=not_f_col,
                    font_color="black",
                )
    pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
    pyvis_graph.write_html(output_name, notebook=False)
    return g


def draw_biomass(
    supermodel: SuperModel,
    output_name: str,
    only_difference=False,
    directed=True,
    n_letter=None,
    wid=None,
    hei=None,
    size=25,
    write_table=True,
    yes_range=1,
    no_range=1,
):
    if not output_name.endswith(".html"):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    color_brewer = get_color_palette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    if write_table:
        output = {
            "Metabolite": [],
            "Type": [],
            "Met_confidence": [],
            "Biomass_confidence": [],
            "Status": [],
        }
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
            if write_table:
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
            if write_table:
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
    pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
    pyvis_graph.write_html(output_name, notebook=False)
    if write_table:
        t_name = re.sub(".html$", ".tsv", output_name)
        output_tb = pd.DataFrame(output)
        output_tb.to_csv(t_name, sep="\t", index=False)
        return g, output_tb
    else:
        return g


def draw_notconv_biomass(
    supermodel: SuperModel,
    output_name: str,
    directed=True,
    n_letter=None,
    wid=None,
    hei=None,
    size=25,
    write_table=True,
    yes_range=1,
    no_range=1,
):
    if not output_name.endswith(".html"):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    g = nx.DiGraph()
    biomass_r = supermodel.reactions.notconverted.get("Biomass")
    reactants = {}
    products = {}
    for typ in supermodel.sources:
        for rea in biomass_r.reactants.get(typ):
            if rea.id not in reactants.keys():
                reactants.update({rea.id: f"{rea.id}\n{typ[:n_letter]}"})
            else:
                if rea.id in rea.annotation.get(typ):
                    reactants.update(
                        {rea.id: f"{reactants.get(rea.id)}{typ[:n_letter]}"}
                    )
        for pro in biomass_r.products.get(typ):
            if pro.id not in products.keys():
                products.update({pro.id: f"{pro.id}\n{typ[:n_letter]}"})
            else:
                products.update({pro.id: f"{products.get(pro.id)}{typ[:n_letter]}"})
    g.add_node(biomass_r.id, shape="box")
    for r in reactants.values():
        g.add_node(r, shape="o")
        g.add_edge(r, biomass_r.id, font_color="black")
    for p in products.values():
        g.add_node(p, shape="o")
        g.add_edge(biomass_r.id, p, font_color="black")

    pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
    pyvis_graph.write_html(output_name, notebook=False)
    return g


def draw_one_synt_path(
    supermodel: SuperModel,
    path: list,
    medium: list,
    met_to_synt: list,
    output_name: str,
    draw_met_not_int=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
    size=25,
):
    if not output_name.endswith(".html"):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    if met_not_int is None:
        met_not_int = {
            "h": 0,
            "h2o": 0,
            "h2": 0,
            "oh1": 0,
            "o2": 0,
            "co2": 0,
            "coa": 0,
            "ppi": 0,
            "pi": 0,
            "amp": 0,
            "adp": 0,
            "atp": 0,
            "cmp": 0,
            "cdp": 0,
            "ctp": 0,
            "gmp": 0,
            "gdp": 0,
            "gtp": 0,
            "ump": 0,
            "udp": 0,
            "utp": 0,
            "nad": 0,
            "nadh": 0,
            "nadp": 0,
            "nadph": 0,
            "dadp": 0,
            "damp": 0,
            "nh3": 0,
            "nh4": 0,
            "fadh2": 0,
            "fad": 0,
            "ac": 0,
            "accoa": 0,
            "h2s": 0,
            "HC00250": 0,
        }
    color_brewer = get_color_palette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    g = nx.DiGraph()
    for r_id in path:
        r = supermodel.reactions.assembly.get(r_id)
        colname_r = define_node_features(
            color_brewer,
            "single_path_r",
            r,
            n_letter,
        )
        g.add_node(
            colname_r[0],
            label=colname_r[1],
            shape="box",
            color=colname_r[2],
            title=colname_r[3],
        )
        for rea in r.reactants.get("assembly"):
            tmp_rea = re.sub("_([cep])$", "", rea.id)
            if (tmp_rea in met_not_int.keys()) and not draw_met_not_int:
                continue
            if (tmp_rea in met_not_int.keys()) and draw_met_not_int:
                colname_rea = define_node_features(
                    color_brewer,
                    "metConnect",
                    rea,
                    n_letter,
                    node_id=f"{tmp_rea} {met_not_int[tmp_rea]}",
                    label=tmp_rea,
                )
                rea_edge = define_edge_features(
                    color_brewer,
                    "single_path_r",
                    "metabolites",
                    rea,
                    f"{tmp_rea} {met_not_int[tmp_rea]}",
                    r,
                    colname_r[1],
                    "reactants",
                    n_letter,
                )
                met_not_int[tmp_rea] = met_not_int[tmp_rea] + 1
            else:
                if tmp_rea in medium:
                    m_pallitra = "metConnect"
                elif tmp_rea in met_to_synt:
                    m_pallitra = "interest"
                else:
                    m_pallitra = "metabolites"
                colname_rea = define_node_features(
                    color_brewer,
                    m_pallitra,
                    rea,
                    n_letter,
                )
                rea_edge = define_edge_features(
                    color_brewer,
                    "single_path_r",
                    "metabolites",
                    rea,
                    colname_rea[1],
                    r,
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
                )
        for pro in r.products.get("assembly"):
            tmp_pro = re.sub("_([cep])$", "", pro.id)
            if (tmp_pro in met_not_int.keys()) and not draw_met_not_int:
                continue
            if (tmp_pro in met_not_int.keys()) and draw_met_not_int:
                colname_pro = define_node_features(
                    color_brewer,
                    "metConnect",
                    pro,
                    n_letter,
                    node_id=f"{tmp_pro} {met_not_int[tmp_pro]}",
                    label=tmp_pro,
                )
                pro_edge = define_edge_features(
                    color_brewer,
                    "single_path_r",
                    "metabolites",
                    pro,
                    f"{tmp_pro} {met_not_int[tmp_pro]}",
                    r,
                    colname_r[1],
                    "products",
                    n_letter,
                )
                met_not_int[tmp_pro] = met_not_int[tmp_pro] + 1
            else:
                if tmp_pro in medium:
                    m_pallitra = "metConnect"
                elif tmp_pro in met_to_synt:
                    m_pallitra = "interest"
                else:
                    m_pallitra = "metabolites"
                colname_pro = define_node_features(
                    color_brewer,
                    m_pallitra,
                    pro,
                    n_letter,
                )
                pro_edge = define_edge_features(
                    color_brewer,
                    "single_path_r",
                    "metabolites",
                    pro,
                    colname_pro[1],
                    r,
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
                )
        if genes:
            g_colname = define_node_features(
                color_brewer,
                "genes",
                r,
                n_letter,
                gene=True,
                and_as_solid=and_as_solid,
            )
            g.add_node(
                g_colname[0],
                shape="box",
                label=g_colname[1],
                color=g_colname[2],
                title=g_colname[3],
            )
            g.add_edge(
                r_id,
                g_colname[0],
                color=g_colname[2],
                font_color="black",
            )
    pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
    pyvis_graph.write_html(output_name, notebook=False)
    return g


def draw_pfba_results(
    path_pfba_out: dict,
    supermodel: SuperModel,
    medium: list,
    write_output_to_folder: str,
    draw_pfba_for_models=None,
    draw_met_not_int=False,
):
    if draw_pfba_for_models is None:
        met_model = {}
        for mod, res in path_pfba_out.items():
            for met in res.keys():
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
            met_model.items(), columns=["Metabolite", "Most confident pfba"]
        ).to_csv(
            f"{write_output_to_folder}/bp_most_confident_pfba_production.tsv",
            index=False,
            sep="\t",
        )
        for m, mod_conf in met_model.items():
            for model_id in mod_conf:
                v = [
                    vv for vv in path_pfba_out[model_id][m] if not vv.startswith("DM_")
                ]
                g = draw_one_synt_path(
                    supermodel,
                    v,
                    medium,
                    [m.removesuffix("_path_pfba")],
                    f"{write_output_to_folder}/{m}_{model_id}.html",
                    draw_met_not_int=draw_met_not_int,
                )
    else:
        for model_id in draw_pfba_for_models:
            for k, v in path_pfba_out[model_id].items():
                v = [vv for vv in v if not vv.startswith("DM_")]
                g = draw_one_synt_path(
                    supermodel,
                    v,
                    medium,
                    [k.removesuffix("_path_pfba")],
                    f"{write_output_to_folder}/{k}_{model_id}.html",
                    draw_met_not_int=draw_met_not_int,
                )


def drawCore(
    supermodel,
    met_not_int,
    colorBrewer,
    name,
    union=False,
    directed=False,
    wid=2000,
    hei=1000,
    Nletter=1,
    union_size=1,
    core_size=None,
):
    g = nx.DiGraph()
    if union:
        attribute = "converted"
        attribute_connect = "union" + str(union_size)
    else:
        if core_size:
            attribute = "core" + str(core_size)
        else:
            attribute = "core" + str(len(supermodel.sources))
        attribute_connect = attribute
    for r in getattr(supermodel.reactions, attribute).values():
        if r.id != "Biomass":
            for rea in r.reactants.get(attribute_connect):
                tmp_rea = (
                    rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                )
                for pro in r.products.get(attribute_connect):
                    tmp_pro = (
                        pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    )
                    if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                        colname_r = define_node_features(
                            colorBrewer, "purples", r.id, supermodel.reactions, Nletter
                        )
                        colname_rea = define_node_features(
                            colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                        )
                        colname_pro = define_node_features(
                            colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                        )
                        cn_rea = define_edge_features(
                            colorBrewer, "purples", rea, r.reactants, Nletter
                        )
                        cn_pro = define_edge_features(
                            colorBrewer, "purples", pro, r.products, Nletter
                        )
                        g.add_node(colname_r[1], shape="box", color=colname_r[0])
                        g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                        g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                        g.add_edge(
                            colname_rea[1],
                            colname_r[1],
                            color=cn_rea[0],
                            font_color="black",
                        )
                        g.add_edge(
                            colname_r[1],
                            colname_pro[1],
                            color=cn_pro[0],
                            font_color="black",
                        )

    pyvis_graph = Network(
        width="{}px".format(wid),
        height="{}px".format(hei),
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawMetSynthesis(
    supermodel,
    met_not_int: list,
    colorBrewer: dict,
    met_of_int: str,
    paths,
    sources: list,
    level: str,
    name_dir: str,
    directed=True,
    wid=2000,
    hei=1000,
    Nletter=1,
):
    g = nx.DiGraph()
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = (
                            rea.id.removesuffix("_c")
                            .removesuffix("_e")
                            .removesuffix("_p")
                        )
                        for pro in r.products.get("core1"):
                            tmp_pro = (
                                pro.id.removesuffix("_c")
                                .removesuffix("_e")
                                .removesuffix("_p")
                            )
                            if (tmp_rea not in met_not_int) & (
                                tmp_pro not in met_not_int
                            ):
                                colname_r = define_node_features(
                                    colorBrewer, "single_path_r", r, Nletter
                                )
                                cn_rea = define_edge_features(
                                    colorBrewer,
                                    "metabolites",
                                    rea,
                                    r,
                                    r.reactants,
                                    Nletter,
                                )
                                cn_pro = define_edge_features(
                                    colorBrewer,
                                    "metabolites",
                                    pro,
                                    r,
                                    r.products,
                                    Nletter,
                                )
                                cn_rea_2 = define_edge_features(
                                    colorBrewer,
                                    "single_path_r",
                                    rea,
                                    r,
                                    r.reactants,
                                    Nletter,
                                )
                                cn_pro_2 = define_edge_features(
                                    colorBrewer,
                                    "single_path_r",
                                    pro,
                                    r,
                                    r.products,
                                    Nletter,
                                )
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = define_node_features(
                                        colorBrewer, "interest", rea, Nletter
                                    )
                                else:
                                    colname_rea = define_node_features(
                                        colorBrewer, "metabolites", rea, Nletter
                                    )
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = define_node_features(
                                        colorBrewer, "interest", pro, Nletter
                                    )
                                else:
                                    colname_pro = define_node_features(
                                        colorBrewer, "metabolites", pro, Nletter
                                    )
                                g.add_node(
                                    colname_r[1], shape="box", color=colname_r[0]
                                )
                                g.add_node(
                                    colname_rea[1], shape="o", color=colname_rea[0]
                                )
                                g.add_node(
                                    colname_pro[1], shape="o", color=colname_pro[0]
                                )
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea_2[0],
                                        font_color="black",
                                    )
                                else:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                    r_done.append(r_id)
        pyvis_graph = Network(
            width="{}px".format(wid),
            height="{}px".format(hei),
            directed=directed,
            notebook=False,
        )
        pyvis_graph.from_nx(g)
        pyvis_graph.show("../Output/" + name_dir + "/" + met_of_int + ".html")
    else:
        print(f"No paths for {met_of_int}, because {paths}")
    return g


def drawMetSynthesis1model(
    g,
    col_r,
    supermodel,
    met_not_int: list,
    colorBrewer: dict,
    met_of_int: str,
    paths,
    sources: list,
    level: str,
    directed=True,
    wid=2000,
    hei=1000,
    Nletter=1,
):
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = (
                            rea.id.removesuffix("_c")
                            .removesuffix("_e")
                            .removesuffix("_p")
                        )
                        for pro in r.products.get("core1"):
                            tmp_pro = (
                                pro.id.removesuffix("_c")
                                .removesuffix("_e")
                                .removesuffix("_p")
                            )
                            if (tmp_rea not in met_not_int) & (
                                tmp_pro not in met_not_int
                            ):
                                colname_r = define_node_features(
                                    colorBrewer, col_r, r, Nletter
                                )
                                cn_rea = define_edge_features(
                                    colorBrewer,
                                    "metabolites",
                                    rea,
                                    r,
                                    r.reactants,
                                    Nletter,
                                )
                                cn_pro = define_edge_features(
                                    colorBrewer,
                                    "metabolites",
                                    pro,
                                    r,
                                    r.products,
                                    Nletter,
                                )
                                cn_rea_2 = define_edge_features(
                                    colorBrewer, col_r, rea, r, r.reactants, Nletter
                                )
                                cn_pro_2 = define_edge_features(
                                    colorBrewer, col_r, pro, r, r.products, Nletter
                                )
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = define_node_features(
                                        colorBrewer, "interest", rea, Nletter
                                    )
                                else:
                                    colname_rea = define_node_features(
                                        colorBrewer, "metabolites", rea, Nletter
                                    )
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = define_node_features(
                                        colorBrewer, "interest", pro, Nletter
                                    )
                                else:
                                    colname_pro = define_node_features(
                                        colorBrewer, "metabolites", pro, Nletter
                                    )
                                g.add_node(
                                    colname_r[1], shape="box", color=colname_r[0]
                                )
                                g.add_node(
                                    colname_rea[1], shape="o", color=colname_rea[0]
                                )
                                g.add_node(
                                    colname_pro[1], shape="o", color=colname_pro[0]
                                )
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea_2[0],
                                        font_color="black",
                                    )
                                else:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                    r_done.append(r_id)
    return g
