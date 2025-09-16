from gemsembler.dbs import (
    get_bigg_network,
    get_kegg_m,
    get_kegg_r,
    get_old_bigg_m,
    get_old_bigg_r,
    get_seed_addit_m,
    get_seed_addit_r,
    get_seed_orig_m,
    get_seed_orig_r,
)


class TestConversionMappings:
    def test_old_bigg(self):
        # Check metabolites
        old_bigg_m = get_old_bigg_m()

        # Check number of mappings
        assert len(old_bigg_m) == 30144

        # Check couple of key-value pairs within the dictionary
        assert old_bigg_m.get("10fthf") == ["10fthf"]
        assert old_bigg_m.get("h2o") == ["h2o"]
        assert old_bigg_m.get("m2mpdol[c]") == ["m2mpdol"]
        assert old_bigg_m.get("lipidA_core_e") == ["lipidA_core_e", "lipidA_core"]

        # Check reactions
        old_bigg_r = get_old_bigg_r()

        # Check number of mappings
        assert len(old_bigg_r) == 32797

        assert old_bigg_r.get("10FTHF5GLUtl") == ["10FTHF5GLUtl"]
        assert old_bigg_r.get("RBK") == ["RBK"]
        assert old_bigg_r.get("DAPAT") == ["DAPAT", "LDAPAT"]

    def test_seed_orig(self):
        # Check metabolites
        seed_orig_m = get_seed_orig_m()

        # Check number of mappings
        assert len(seed_orig_m) == 2652

        # Check couble of key-value pairs within the dictionary
        assert seed_orig_m.get("cpd00001") == ["h2o", "oh1"]
        assert seed_orig_m.get("cpd00200") == ["4mop"]
        assert seed_orig_m.get("cpd00110") == ["focytC", "focytc", "focytcc553"]

        # Check reactions
        seed_orig_r = get_seed_orig_r()

        # Check number of mappings
        assert len(seed_orig_r) == 7331

        # Check couble of key-value pairs within the dictionary
        assert seed_orig_r.get("rxn00001") == ["IPP1", "PPA", "PPA_1", "PPAm"]
        assert seed_orig_r.get("rxn00248") == ["MDH", "MDH1", "MDHi2", "MDHm", "MDHp"]
        assert seed_orig_r.get("rxn00456") == ["METGL"]

    def test_seed_addit(self):
        # Check metabolites
        seed_addit_m = get_seed_addit_m()

        # Check number of mappings
        assert len(seed_addit_m) == 2910

        # Check couble of key-value pairs within the dictionary
        assert seed_addit_m.get("cpd15275") == ["oh1"]
        assert seed_addit_m.get("cpd23593") == ["cdpdhdecg"]

        # Check reactions
        seed_addit_r = get_seed_addit_r()

        # Check number of mappings
        assert len(seed_addit_r) == 7741

        # Check couble of key-value pairs within the dictionary
        assert seed_addit_r.get("rxn22163") == [
            "CRBNTD",
            "DNADRAIN",
            "H2CO3D2",
            "H2CO3D2m",
            "HMR_5409",
            "NH3c",
            "NH4DISg",
            "RE2594C",
        ]
        assert seed_addit_r.get("rxn35210") == ["GCLDH"]
        assert seed_addit_r.get("rxn00184") == [
            "GDH_nadp",
            "GLUDy",
            "GLUDy_1",
            "GLUDym",
        ]

    def test_kegg(self):
        # Check metabolites
        kegg_m = get_kegg_m()

        # Check number of mappings
        assert len(kegg_m) == 2160

        # Check couble of key-value pairs within the dictionary
        assert kegg_m.get("C01328") == ["oh1"]
        assert kegg_m.get("C15547") == ["14dh2napcoa", "14dhncoa", "dhncoa"]

        # Check reactions
        kegg_r = get_kegg_r()

        # Check number of mappings
        assert len(kegg_r) == 1739

        # Check couble of key-value pairs within the dictionary
        assert kegg_r.get("R00253") == ["GALh", "GALm", "GLNS", "GLNS_1"]
        assert kegg_r.get("R03057") == ["LTA4H"]
        assert kegg_r.get("R00703") == ["LDH_L", "LDH_Lm", "LLDHh", "r0173"]

    def test_bigg_network(self):
        bigg_net = get_bigg_network()

        # Check number of reactions
        assert len(bigg_net) == 27946

        # Check whether specific key exists in the bigg_net dict
        assert "10fthf5glu_c adp_c h2o_c h_c pi_c<->10fthf_c atp_c glu__L_c" in bigg_net

        # Check whether key-value pair matches within bigg_net dict
        assert (
            bigg_net["10fthf5glu_c adp_c h2o_c h_c pi_c<->10fthf_c atp_c glu__L_c"]
            == "FPGS7_1"
        )
