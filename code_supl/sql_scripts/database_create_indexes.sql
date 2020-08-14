----Training sets
CREATE INDEX ON chemical_interactions_v2.predictor_1 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1 (ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.predictor_2 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2 (ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.predictor_1_high_scored (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_high_scored (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_high_scored (ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.predictor_2_high_scored (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_high_scored (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_high_scored (ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.predictor_1_independent_validation (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_independent_validation (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_independent_validation (ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.predictor_2_independent_validation (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_independent_validation (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_independent_validation (ensp_1, ensp_2);


----Drugs
CREATE INDEX ON chemical_interactions_v2.all_tested_drugs (stitch_pubchem_cid);

CREATE INDEX ON chemical_interactions_v2.drug_atc (cid_1);
CREATE INDEX ON chemical_interactions_v2.drug_atc (cid_2);
CREATE INDEX ON chemical_interactions_v2.drug_atc (cid_1, cid_2);

CREATE INDEX ON chemical_interactions_v2.drug_side_effect (cid_1);
CREATE INDEX ON chemical_interactions_v2.drug_side_effect (cid_2);
CREATE INDEX ON chemical_interactions_v2.drug_side_effect (cid_1, cid_2);


----Targets
CREATE INDEX ON chemical_interactions_v2.target_alias(ensp);
CREATE INDEX ON chemical_interactions_v2.target_alias(uniprot_id);
CREATE INDEX ON chemical_interactions_v2.target_alias(entrez_kegg_id);


----FxnAss
CREATE INDEX ON chemical_interactions_v2.string_topo(ensp_1);
CREATE INDEX ON chemical_interactions_v2.string_topo(ensp_2);
CREATE INDEX ON chemical_interactions_v2.string_topo(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.string_topo_eb(ensp_1);
CREATE INDEX ON chemical_interactions_v2.string_topo_eb(ensp_2);
CREATE INDEX ON chemical_interactions_v2.string_topo_eb(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.string_topo_nsp(ensp_1);
CREATE INDEX ON chemical_interactions_v2.string_topo_nsp(ensp_2);
CREATE INDEX ON chemical_interactions_v2.string_topo_nsp(ensp_1, ensp_2);


----GeneEss
CREATE INDEX ON chemical_interactions_v2.gene_essentiality(ensp);


----GeneExp
CREATE INDEX ON chemical_interactions_v2.gene_coexpression(ensp_1);
CREATE INDEX ON chemical_interactions_v2.gene_coexpression(ensp_2);
CREATE INDEX ON chemical_interactions_v2.gene_coexpression(ensp_1, ensp_2);


----GetInt
CREATE INDEX ON chemical_interactions_v2.getint_topo(ensp_1);
CREATE INDEX ON chemical_interactions_v2.getint_topo(ensp_2);
CREATE INDEX ON chemical_interactions_v2.getint_topo(ensp_1, ensp_2);


CREATE INDEX ON chemical_interactions_v2.getint_topo_eb(ensp_1);
CREATE INDEX ON chemical_interactions_v2.getint_topo_eb(ensp_2);
CREATE INDEX ON chemical_interactions_v2.getint_topo_eb(ensp_1, ensp_2);


CREATE INDEX ON chemical_interactions_v2.getint_topo_nsp(ensp_1);
CREATE INDEX ON chemical_interactions_v2.getint_topo_nsp(ensp_2);
CREATE INDEX ON chemical_interactions_v2.getint_topo_nsp(ensp_1, ensp_2);


----PPI
CREATE INDEX ON chemical_interactions_v2.biogrid_topo(ensp_1);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo(ensp_2);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo(ensp_1, ensp_2);


CREATE INDEX ON chemical_interactions_v2.biogrid_topo_eb(ensp_1);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo_eb(ensp_2);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo_eb(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.biogrid_topo_nsp(ensp_1);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo_nsp(ensp_2);
CREATE INDEX ON chemical_interactions_v2.biogrid_topo_nsp(ensp_1, ensp_2);


----GO
CREATE INDEX ON chemical_interactions_v2.go_all(ensp_1);
CREATE INDEX ON chemical_interactions_v2.go_all(ensp_2);
CREATE INDEX ON chemical_interactions_v2.go_all(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.go_bp(ensp_1);
CREATE INDEX ON chemical_interactions_v2.go_bp(ensp_2);
CREATE INDEX ON chemical_interactions_v2.go_bp(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.go_cc(ensp_1);
CREATE INDEX ON chemical_interactions_v2.go_cc(ensp_2);
CREATE INDEX ON chemical_interactions_v2.go_cc(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.go_mf(ensp_1);
CREATE INDEX ON chemical_interactions_v2.go_mf(ensp_2);
CREATE INDEX ON chemical_interactions_v2.go_mf(ensp_1, ensp_2);

CREATE INDEX ON chemical_interactions_v2.go_slim(ensp_1);
CREATE INDEX ON chemical_interactions_v2.go_slim(ensp_2);
CREATE INDEX ON chemical_interactions_v2.go_slim(ensp_1, ensp_2);


---- Phylo
CREATE INDEX ON chemical_interactions_v2.phylo(ensp_1);
CREATE INDEX ON chemical_interactions_v2.phylo(ensp_2);
CREATE INDEX ON chemical_interactions_v2.phylo(ensp_1, ensp_2);





