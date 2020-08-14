CREATE TABLE chemical_interactions_v2.predictor_all_try2 (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_all_try2 (
SELECT pcl_1.cid cid_1, pcl_2.cid cid_2, ensp_1, ensp_2
FROM chemical_interactions_v2.biogrid_topo
JOIN chemical_interactions_v2.biogrid_topo_eb USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.gene_coexpression USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.gene_essentiality gene_ess_1 ON (gene_ess_1.ensp = ensp_1)
JOIN chemical_interactions_v2.gene_essentiality gene_ess_2 ON (gene_ess_2.ensp = ensp_2)
JOIN chemical_interactions_v2.getint_topo USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.getint_topo_eb USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.go_all USING (ensp_1, ensp_2)
LEFT JOIN chemical_interactions_v2.go_bp USING (ensp_1, ensp_2)
LEFT JOIN chemical_interactions_v2.go_cc USING (ensp_1, ensp_2)
LEFT JOIN chemical_interactions_v2.go_mf USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.phylo USING (ensp_1, ensp_2)
JOIN chemical_interactions_v2.string_topo USING (ensp_1, ensp_2)
LEFT JOIN chemical_interactions_v2.string_topo_eb USING (ensp_1, ensp_2)
JOIN stitch.protein_chemical_links_human_nostereo_hc pcl_1 ON (pcl_1.ensp = ensp_1)
JOIN stitch.protein_chemical_links_human_nostereo_hc pcl_2 ON (pcl_2.ensp = ensp_2));

CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_all_try2 (ensp_1, ensp_2);




CREATE TABLE chemical_interactions_v2.predictor_2_all_try2 (
    ensp_1 int,
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_2_all_try2 (
SELECT DISTINCT ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_all)

CREATE INDEX ON chemical_interactions_v2.predictor_2_all_try2 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_try2 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_try2 (ensp_1, ensp_2);




CREATE TABLE chemical_interactions_v2.predictor_1_all_try2 (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_1_all_try2 (
SELECT DISTINCT cid_1, cid_2, ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_all
JOIN chemical_interactions_v2.drug_atc USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_chemical_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_side_effect USING (cid_1, cid_2) );

CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_try2 (ensp_1, ensp_2);





