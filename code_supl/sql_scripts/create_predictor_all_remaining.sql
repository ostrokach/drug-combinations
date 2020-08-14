DISCARD TEMP;

CREATE TEMPORARY TABLE cids_to_exclude_1 (cid int);

INSERT INTO cids_to_exclude_1
SELECT DISTINCT cid_1 cid
FROM chemical_interactions_v2.predictor_1
UNION
SELECT DISTINCT cid_2 cid
FROM chemical_interactions_v2.predictor_1
UNION
SELECT DISTINCT cid_1 cid
FROM chemical_interactions_v2.predictor_1_independent_validation
UNION
SELECT DISTINCT cid_2 cid
FROM chemical_interactions_v2.predictor_1_independent_validation;

CREATE INDEX cids_to_exlude_idx ON cids_to_exclude_1 (cid);




CREATE TEMPORARY TABLE ensps_to_exclude_1 (ensp int);

INSERT INTO ensps_to_exclude_1
SELECT DISTINCT ensp_1 ensp
FROM chemical_interactions_v2.predictor_1
UNION
SELECT DISTINCT ensp_2 ensp
FROM chemical_interactions_v2.predictor_1
UNION
SELECT DISTINCT ensp_1 ensp
FROM chemical_interactions_v2.predictor_1_independent_validation
UNION
SELECT DISTINCT ensp_2 ensp
FROM chemical_interactions_v2.predictor_1_independent_validation;

CREATE INDEX ensps_to_exclude_idx ON ensps_to_exclude_1 (ensp);




CREATE TEMPORARY TABLE ensps_to_exclude_2 (ensp int);

INSERT INTO ensps_to_exclude_2
SELECT DISTINCT ensp_1 ensp
FROM chemical_interactions_v2.predictor_2
UNION
SELECT DISTINCT ensp_2 ensp
FROM chemical_interactions_v2.predictor_2
UNION
SELECT DISTINCT ensp_1 ensp
FROM chemical_interactions_v2.predictor_2_independent_validation
UNION
SELECT DISTINCT ensp_2 ensp
FROM chemical_interactions_v2.predictor_2_independent_validation;

CREATE INDEX ensps_to_exclude_idx ON ensps_to_exclude_2 (ensp);




DROP TABLE IF EXISTS chemical_interactions_v2.predictor_1_all_unused;

CREATE TABLE chemical_interactions_v2.predictor_1_all_unused (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);
    
INSERT INTO chemical_interactions_v2.predictor_1_all_unused (
SELECT cid_1, cid_2, ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_1_all
WHERE NOT EXISTS (SELECT 1 FROM cids_to_exclude_1 WHERE cid=cid_1 OR cid=cid_2)
AND NOT EXISTS (SELECT 1 FROM ensps_to_exclude_1 WHERE ensp=ensp_1 OR ensp=ensp_2));

CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (cid_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (cid_1, cid_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_1_all_unused (ensp_1, ensp_2);




DROP TABLE IF EXISTS chemical_interactions_v2.predictor_2_all_unused;

CREATE TABLE chemical_interactions_v2.predictor_2_all_unused (
    ensp_1 int,
    ensp_2 int);
    
INSERT INTO chemical_interactions_v2.predictor_2_all_unused (
SELECT ensp_1, ensp_2
FROM chemical_interactions_v2.predictor_2_all
WHERE NOT EXISTS (SELECT 1 FROM ensps_to_exclude_2 WHERE ensp=ensp_1 OR ensp=ensp_2));

CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_2_all_unused (ensp_1, ensp_2);



