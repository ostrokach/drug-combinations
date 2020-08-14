CREATE TABLE chemical_interactions_v2.predictor_2_all (
    ensp_1 int, 
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_2_all (
SELECT DISTINCT a.ensp_1, a.ensp_2
FROM chemical_interactions_v2.biogrid_topo a
WHERE EXISTS (SELECT 1 FROM chemical_interactions_v2.biogrid_topo_eb WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.gene_coexpression WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.getint_topo WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.getint_topo_eb WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.go_all WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.phylo WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.string_topo WHERE ensp_1=a.ensp_1 AND ensp_2=a.ensp_2)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.gene_essentiality ge_1 WHERE ge_1.ensp=a.ensp_1)
AND EXISTS (SELECT 1 FROM chemical_interactions_v2.gene_essentiality ge_2 WHERE ge_2.ensp=a.ensp_2)
AND EXISTS (SELECT 1 FROM stitch.protein_chemical_links_human_nostereo_hc WHERE ensp=a.ensp_1)
AND EXISTS (SELECT 1 FROM stitch.protein_chemical_links_human_nostereo_hc WHERE ensp=a.ensp_2) );

CREATE INDEX ON chemical_interactions_v2.predictor_2_all (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all (ensp_1, ensp_2);



DROP TABLE chemical_interactions_v2.predictor_1_all;

CREATE TABLE chemical_interactions_v2.predictor_1_all (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_1_all (
SELECT DISTINCT pcl_1.cid, pcl_2.cid, pcl_1.ensp, pcl_2.ensp
FROM 
(SELECT cid_1, cid_2
FROM chemical_interactions_v2.drug_atc_similarity a
WHERE EXISTS (SELECT 1 FROM chemical_interactions_v2.drug_side_effect_similarity WHERE cid_1=a.cid_1 AND cid_2=a.cid_2)) a
--AND EXISTS (SELECT 1 FROM chemical_interactions_v2.drug_chemical_similarity WHERE cid_1=a.cid_1 AND cid_2=a.cid_2)
JOIN stitch.protein_chemical_links_human_nostereo_hc pcl_1 ON (pcl_1.cid=a.cid_1)
JOIN stitch.protein_chemical_links_human_nostereo_hc pcl_2 ON (pcl_2.cid=a.cid_2)
WHERE EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_2_all WHERE ensp_1=pcl_1.ensp AND ensp_2=pcl_2.ensp));

CREATE INDEX ON chemical_interactions_v2.predictor_1_all (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all (ensp_1, ensp_2);

