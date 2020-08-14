-- Predictor 1
CREATE TABLE chemical_interactions_v2.predictor_1_all_unused_pairs_3 AS
SELECT cid_1, cid_2, ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_1_all a
WHERE 
NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1 b WHERE b.cid_1=a.cid_1 AND b.cid_2=a.cid_2) AND
NOT EXISTS (SELECT 1 FROM chemical_interactions_v2.predictor_1_independent_validation b WHERE b.cid_1=a.cid_1 AND b.cid_2=a.cid_2);

CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused_pairs_3 (ensp_1, ensp_2);



-- Predictor 2
CREATE TABLE chemical_interactions_v2.predictor_2_all_unused_pairs_3 AS
SELECT ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_2_all a
WHERE 
(ensp_1, ensp_2) not in (select ensp_1, ensp_2 from chemical_interactions_v2.predictor_2) and
(ensp_1, ensp_2) not in (select ensp_1, ensp_2 from chemical_interactions_v2.predictor_2_independent_validation);

CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs_3 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs_3 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused_pairs_3 (ensp_1, ensp_2);
