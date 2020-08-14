DROP TABLE IF EXISTS chemical_interactions_v2.predictor_1_all_unused_pairs;

CREATE TABLE chemical_interactions_v2.predictor_1_all_unused_pairs (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);
    
INSERT INTO chemical_interactions_v2.predictor_1_all_unused_pairs (
SELECT cid_1, cid_2, ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_1_all a
WHERE NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1 b WHERE b.cid_1=a.cid_1 OR b.cid_2=a.cid_2)
AND NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1_independent_validation b WHERE b.cid_1=a.cid_1 OR b.cid_2=a.cid_2)
AND NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1 b WHERE b.ensp_1=a.ensp_1 OR b.ensp_2=a.ensp_2)
AND NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1_independent_validation b WHERE b.ensp_1=a.ensp_1 OR b.ensp_2=a.ensp_2));

CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs (ensp_1, ensp_2);




DROP TABLE IF EXISTS chemical_interactions_v2.predictor_2_all_unused_pairs;

CREATE TABLE chemical_interactions_v2.predictor_2_all_unused_pairs (
    ensp_1 int,
    ensp_2 int);
    
INSERT INTO chemical_interactions_v2.predictor_2_all_unused_pairs (
SELECT ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_2_all a
WHERE NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_2 b WHERE b.ensp_1=a.ensp_1 OR b.ensp_2=a.ensp_2)
AND NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_2_independent_validation b WHERE b.ensp_1=a.ensp_1 OR b.ensp_2=a.ensp_2));

CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs (ensp_1, ensp_2);



