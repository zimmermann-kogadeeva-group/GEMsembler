from importlib.resources import files

from gemsembler import load_sbml_model
from gemsembler.conversion import (
    ConvAgora,
    ConvBase,
    ConvCarveme,
    ConvGapseq,
    ConvModelseed,
)
from gemsembler.data import BU


class TestConversion:
    def test_agora(self):
        # Open model and get a metabolite and a reaction to check
        model = load_sbml_model(files(BU) / "BU_agora.xml.gz")

        # Create a conversion object
        conv = ConvAgora()

        # Convert a metabolite
        mb = model.metabolites.get_by_id("2dmmq8[c]")
        conv_mb = conv.convert_metabolite(mb)

        # Check the converion of the metabolite
        assert conv_mb.compartments == ["c"]
        assert conv_mb.highest == ["2dmmq8_c"]
        assert conv_mb.level == "annot_and_main"
        assert conv_mb.annot_and_main == ["2dmmq8_c"]
        assert conv_mb.annot == []
        assert conv_mb.main == []
        assert conv_mb.addit == []
        assert conv_mb.pattern == []
        assert conv_mb.no_conv == []

        # Convert a reaction
        reac = model.reactions.get_by_id("ADNCNT3tc")
        conv_reac = conv.convert_reaction(reac)

        # Check the converion of the reaction
        assert set(conv_reac.compartments) == {"c", "e"}
        assert conv_reac.highest == ["ADNt2"]
        assert conv_reac.level == "main"
        assert conv_reac.annot_and_main == []
        assert conv_reac.annot == []
        assert conv_reac.main == ["ADNt2"]
        assert conv_reac.addit == []
        assert conv_reac.pattern == []
        assert conv_reac.no_conv == []

        # Convert the whole model and count number of all conversions
        conv_mbs, conv_reacs = conv.convert_model(model).values()

        # Metabolites
        assert len(conv_mbs) == 1500

        num_converted = len([x for x in conv_mbs.values() if x.highest])
        assert num_converted == 892

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_mbs.values()])
        assert num_annot_and_main == 799

        num_annot = sum([len(x.annot) for x in conv_mbs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_mbs.values()])
        assert num_main == 5

        num_addit = sum([len(x.addit) for x in conv_mbs.values()])
        assert num_addit == 109

        num_pattern = sum([len(x.pattern) for x in conv_mbs.values()])
        assert num_pattern == 21

        num_no_conv = sum([len(x.no_conv) for x in conv_mbs.values()])
        assert num_no_conv == 0

        # Reactions
        assert len(conv_reacs) == 2418

        num_converted = len([x for x in conv_reacs.values() if x.highest])
        assert num_converted == 849

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_reacs.values()])
        assert num_annot_and_main == 713

        num_annot = sum([len(x.annot) for x in conv_reacs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_reacs.values()])
        assert num_main == 377

        num_addit = sum([len(x.addit) for x in conv_reacs.values()])
        assert num_addit == 245

        num_pattern = sum([len(x.pattern) for x in conv_reacs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_reacs.values()])
        assert num_no_conv == 0

    def test_gapseq(self):
        # Open model and get a metabolite and a reaction to check
        model = load_sbml_model(files(BU) / "BU_gapseq.xml.gz")

        # Create a conversion object
        conv = ConvGapseq()

        # Convert a metabolite
        mb = model.metabolites.get_by_id("cpd00843_c0")
        conv_mb = conv.convert_metabolite(mb)

        # Check the converion of the metabolite
        assert conv_mb.compartments == ["c"]
        assert conv_mb.highest == ["2h3oppan_c"]
        assert conv_mb.level == "annot_and_main"
        assert conv_mb.annot_and_main == ["2h3oppan_c"]
        assert set(conv_mb.annot) == set(["2h3opp_c", "hop_c"])
        assert conv_mb.main == []
        assert conv_mb.addit == []
        assert conv_mb.pattern == []
        assert conv_mb.no_conv == []

        # Convert a reaction
        reac = model.reactions.get_by_id("rxn01301_c0")
        conv_reac = conv.convert_reaction(reac)

        # Check the converion of the reaction
        assert conv_reac.compartments == ["c"]
        assert conv_reac.highest == ["HSDxi"]
        assert conv_reac.level == "annot_and_main"
        assert conv_reac.annot_and_main == ["HSDxi"]
        assert conv_reac.annot == []
        assert conv_reac.main == []
        assert set(conv_reac.addit) == {"HSDH_h", "HSDH_m"}
        assert conv_reac.pattern == []
        assert conv_reac.no_conv == []

        # Convert the whole model and count number of all conversions
        conv_mbs, conv_reacs = conv.convert_model(model).values()

        # Metabolites
        assert len(conv_mbs) == 1520

        num_converted = len([x for x in conv_mbs.values() if x.highest])
        assert num_converted == 1219

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_mbs.values()])
        assert num_annot_and_main == 1002

        num_annot = sum([len(x.annot) for x in conv_mbs.values()])
        assert num_annot == 114

        num_main = sum([len(x.main) for x in conv_mbs.values()])
        assert num_main == 179

        num_addit = sum([len(x.addit) for x in conv_mbs.values()])
        assert num_addit == 94

        num_pattern = sum([len(x.pattern) for x in conv_mbs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_mbs.values()])
        assert num_no_conv == 0

        # Reactions
        assert len(conv_reacs) == 1891

        num_converted = len([x for x in conv_reacs.values() if x.highest])
        assert num_converted == 1213

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_reacs.values()])
        assert num_annot_and_main == 828

        num_annot = sum([len(x.annot) for x in conv_reacs.values()])
        assert num_annot == 631

        num_main = sum([len(x.main) for x in conv_reacs.values()])
        assert num_main == 509

        num_addit = sum([len(x.addit) for x in conv_reacs.values()])
        assert num_addit == 702

        num_pattern = sum([len(x.pattern) for x in conv_reacs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_reacs.values()])
        assert num_no_conv == 0

    def test_modelseed(self):
        # Open model and get a metabolite and a reaction to check
        model = load_sbml_model(files(BU) / "BU_modelSEED.sbml.gz")

        # Create a conversion object
        conv = ConvModelseed()

        # Convert a metabolite
        mb = model.metabolites.get_by_id("cpd00013_e0")
        conv_mb = conv.convert_metabolite(mb)

        # Check the converion of the metabolite
        assert conv_mb.compartments == ["e"]
        assert set(conv_mb.highest) == set(["nh4_e", "nh3_e"])
        assert conv_mb.level == "main"
        assert conv_mb.annot_and_main == []
        assert conv_mb.annot == []
        assert set(conv_mb.main) == set(["nh4_e", "nh3_e"])
        assert conv_mb.addit == []
        assert conv_mb.pattern == []
        assert conv_mb.no_conv == []

        # Convert a reaction
        reac = model.reactions.get_by_id("rxn03108_c0")
        conv_reac = conv.convert_reaction(reac)

        # Check the converion of the reaction
        assert conv_reac.compartments == ["c"]
        assert conv_reac.highest == ["PMPK"]
        assert conv_reac.level == "main"
        assert conv_reac.annot_and_main == []
        assert conv_reac.annot == []
        assert conv_reac.main == ["PMPK"]
        assert conv_reac.addit == ["PMPK_h"]
        assert conv_reac.pattern == []
        assert conv_reac.no_conv == []

        # Convert the whole model and count number of all conversions
        conv_mbs, conv_reacs = conv.convert_model(model).values()

        # Metabolites
        assert len(conv_mbs) == 1332

        num_converted = len([x for x in conv_mbs.values() if x.highest])
        assert num_converted == 969

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_mbs.values()])
        assert num_annot_and_main == 0

        num_annot = sum([len(x.annot) for x in conv_mbs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_mbs.values()])
        assert num_main == 966

        num_addit = sum([len(x.addit) for x in conv_mbs.values()])
        assert num_addit == 109

        num_pattern = sum([len(x.pattern) for x in conv_mbs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_mbs.values()])
        assert num_no_conv == 0

        # Reactions
        assert len(conv_reacs) == 1222

        num_converted = len([x for x in conv_reacs.values() if x.highest])
        assert num_converted == 722

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_reacs.values()])
        assert num_annot_and_main == 0

        num_annot = sum([len(x.annot) for x in conv_reacs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_reacs.values()])
        assert num_main == 869

        num_addit = sum([len(x.addit) for x in conv_reacs.values()])
        assert num_addit == 654

        num_pattern = sum([len(x.pattern) for x in conv_reacs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_reacs.values()])
        assert num_no_conv == 0

    def test_carveme(self):
        # Open model and get a metabolite and a reaction to check
        model = load_sbml_model(files(BU) / "BU_carveme_hom.xml.gz")

        # Create a conversion object
        conv = ConvCarveme()

        # Convert a metabolite
        mb = model.metabolites.get_by_id("dhlam_c")
        conv_mb = conv.convert_metabolite(mb)

        # Check the converion of the metabolite
        assert conv_mb.compartments == ["c"]
        assert conv_mb.highest == ["dhlam_c"]
        assert conv_mb.level == "main"
        assert conv_mb.annot_and_main == []
        assert conv_mb.annot == []
        assert conv_mb.main == ["dhlam_c"]
        assert conv_mb.addit == []
        assert conv_mb.pattern == []
        assert conv_mb.no_conv == []

        # Convert a reaction
        reac = model.reactions.get_by_id("AACPS3")
        conv_reac = conv.convert_reaction(reac)

        # Check the converion of the reaction
        assert conv_reac.compartments == ["c"]
        assert conv_reac.highest == ["AACPS3"]
        assert conv_reac.level == "main"
        assert conv_reac.annot_and_main == []
        assert conv_reac.annot == []
        assert conv_reac.main == ["AACPS3"]
        assert conv_reac.addit == []
        assert conv_reac.pattern == []
        assert conv_reac.no_conv == []

        # Convert the whole model and count number of all conversions
        conv_mbs, conv_reacs = conv.convert_model(model).values()

        # Metabolites
        assert len(conv_mbs) == 1244

        num_converted = len([x for x in conv_mbs.values() if x.highest])
        assert num_converted == 1244

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_mbs.values()])
        assert num_annot_and_main == 0

        num_annot = sum([len(x.annot) for x in conv_mbs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_mbs.values()])
        assert num_main == 1244

        num_addit = sum([len(x.addit) for x in conv_mbs.values()])
        assert num_addit == 0

        num_pattern = sum([len(x.pattern) for x in conv_mbs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_mbs.values()])
        assert num_no_conv == 0

        # Reactions
        assert len(conv_reacs) == 1849

        num_converted = len([x for x in conv_reacs.values() if x.highest])
        assert num_converted == 1839

        num_annot_and_main = sum([len(x.annot_and_main) for x in conv_reacs.values()])
        assert num_annot_and_main == 0

        num_annot = sum([len(x.annot) for x in conv_reacs.values()])
        assert num_annot == 0

        num_main = sum([len(x.main) for x in conv_reacs.values()])
        assert num_main == 1839

        num_addit = sum([len(x.addit) for x in conv_reacs.values()])
        assert num_addit == 0

        num_pattern = sum([len(x.pattern) for x in conv_reacs.values()])
        assert num_pattern == 0

        num_no_conv = sum([len(x.no_conv) for x in conv_reacs.values()])
        assert num_no_conv == 0
