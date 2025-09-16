from gemsembler import GatheredModels, bu_example
from gemsembler.selection import run_selection


class TestSelection:
    def test_first_stage(self):
        g = GatheredModels()
        for model in bu_example:
            g.add_model(**model)

        sel_mbs = run_selection(g.same_db_models, g.converted_metabolites, "highest")

        # Check the number of metabolites is consistant from `run_selection`
        assert len(sel_mbs["carveme_BU"]) == 1244
        assert len(sel_mbs["agora_BU"]) == 1500
        assert len(sel_mbs["gapseq_BU"]) == 1520
        assert len(sel_mbs["modelseed_BU"]) == 1271

        # Check individual metabolite from carveme model
        sel_mb = sel_mbs["carveme_BU"]["10fthf_c"]
        assert sel_mb.compartments == ["c"]
        assert sel_mb.consistent == "Yes"
        assert sel_mb.highest_consistent == ["10fthf_c"]
        assert sel_mb.from_one_id is True
        assert sel_mb.to_one_id is True
        assert sel_mb.from_many_other_ids == []
        assert sel_mb.in_other_models == {}

        # Check individual metabolite from agora model
        sel_mb = sel_mbs["agora_BU"]["gam1p[c]"]
        assert sel_mb.compartments == ["c"]
        assert sel_mb.consistent == "Yes"
        assert sel_mb.highest_consistent == ["gam1p_c"]
        assert sel_mb.from_one_id is True
        assert sel_mb.to_one_id is True
        assert sel_mb.from_many_other_ids == []
        assert sel_mb.in_other_models == {}

        # Check individual metabolite from modelseed model
        sel_mb = sel_mbs["modelseed_BU"]["cpd02902_c0"]
        assert sel_mb.compartments == ["c"]
        assert sel_mb.consistent == "Not converted"
        assert sel_mb.highest_consistent == []
        assert sel_mb.from_one_id is False
        assert sel_mb.to_one_id is None
        assert len(sel_mb.from_many_other_ids) == 355
        assert sel_mb.in_other_models == {}

        # Check individual metabolite from gapseq model
        sel_mb = sel_mbs["gapseq_BU"]["cpd11477_c0"]
        assert sel_mb.compartments == ["c"]
        assert sel_mb.consistent == "Yes"
        assert sel_mb.highest_consistent == ["tpalm2eACP_c"]
        assert sel_mb.from_one_id is True
        assert sel_mb.to_one_id is True
        assert sel_mb.from_many_other_ids == []
        assert sel_mb.in_other_models == {"modelseed_BU": [True, True]}
