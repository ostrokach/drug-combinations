CREATE TEMPORARY TABLE cids_to_include (cid int);

WITH temp AS (
    SELECT DISTINCT cid_1, cid_2
    FROM chemical_interactions_v2.drug_atc
    JOIN chemical_interactions_v2.drug_chemical_similarity USING (cid_1, cid_2)
    JOIN chemical_interactions_v2.drug_side_effect USING (cid_1, cid_2) )
INSERT INTO cids_to_include (
    SELECT DISTINCT cid_1 cid FROM temp
    UNION
    SELECT DISTINCT cid_2 cid FROM temp);

CREATE INDEX ON cids_to_include (cid);




CREATE TEMPORARY TABLE ensps_to_include (ensp int);

INSERT INTO ensps_to_include (
    SELECT DISTINCT ensp_1 ensp 
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
    UNION
    SELECT DISTINCT ensp_2 ensp 
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
    );

CREATE INDEX ON ensps_to_include (ensp);




CREATE TEMPORARY TABLE protein_chemical_links_to_include_2 (
cid integer,
ensp integer);

INSERT INTO protein_chemical_links_to_include_2 (
SELECT DISTINCT cid, ensp
FROM stitch.protein_chemical_links_human_nostereo
JOIN ensps_to_include USING (ensp)
WHERE rank < 3);

CREATE INDEX ON protein_chemical_links_to_include_2 (cid);
CREATE INDEX ON protein_chemical_links_to_include_2 (ensp);
CREATE INDEX ON protein_chemical_links_to_include_2 (cid, ensp);




CREATE TEMPORARY TABLE protein_chemical_links_to_include_2ensponly (ensp integer);

INSERT INTO protein_chemical_links_to_include_2ensponly (
SELECT DISTINCT ensp
FROM protein_chemical_links_to_include_2);

CREATE INDEX ON protein_chemical_links_to_include_2ensponly (ensp);



DROP TABLE IF EXISTS chemical_interactions_v2.predictor_all_2_v2;

CREATE TABLE chemical_interactions_v2.predictor_all_2_v2 (
    ensp_1 int,
    ensp_2 int);

INSERT INTO chemical_interactions_v2.predictor_all_2_v2 (
    (SELECT
        p_1.ensp ensp_1, 
        p_2.ensp ensp_2
    FROM protein_chemical_links_to_include_2ensponly p_1
    JOIN protein_chemical_links_to_include_2ensponly p_2 ON (p_1.ensp < p_2.ensp))
    UNION
    (SELECT
        p_2.ensp ensp_1, 
        p_1.ensp ensp_2
    FROM protein_chemical_links_to_include_2ensponly p_1
    JOIN protein_chemical_links_to_include_2ensponly p_2 ON (p_1.ensp > p_2.ensp)));

CREATE INDEX ON chemical_interactions_v2.predictor_all_2_v2 (ensp_1);
CREATE INDEX ON chemical_interactions_v2.predictor_all_2_v2 (ensp_2);
CREATE INDEX ON chemical_interactions_v2.predictor_all_2_v2 (ensp_1, ensp_2);




CREATE TEMPORARY TABLE protein_chemical_links_to_include_1 (
cid integer,
ensp integer);

INSERT INTO protein_chemical_links_to_include_1 (
SELECT DISTINCT cid, ensp
FROM stitch.protein_chemical_links_to_include_2
JOIN cids_to_include USING (cid) );

CREATE INDEX ON protein_chemical_links_to_include_1 (cid);
CREATE INDEX ON protein_chemical_links_to_include_1 (ensp);
CREATE INDEX ON protein_chemical_links_to_include_1 (cid, ensp);




DROP TABLE IF EXISTS chemical_interactions_v2.predictor_all_1;

CREATE TABLE chemical_interactions_v2.predictor_all_1 (
    cid_1 int,
    cid_2 int,
    ensp_1 int,
    ensp_2 int);

WITH p_pair AS (
    SELECT DISTINCT
        p_1.cid cid_1, 
        p_2.cid cid_2, 
        p_1.ensp ensp_1, 
        p_2.ensp ensp_2
    FROM protein_chemical_links_to_include_1 p_1
    JOIN protein_chemical_links_to_include_1 p_2 ON (p_1.cid < p_2.cid and p_1.ensp <> p_2.ensp) )
INSERT INTO chemical_interactions_v2.predictor_all_1 (
    (SELECT DISTINCT
        p_pair.cid_1 cid_1, 
        p_pair.cid_2 cid_2, 
        p_pair.ensp_1 ensp_1, 
        p_pair.ensp_2 ensp_2
    FROM p_pair
    WHERE ensp_1 < ensp_2)
    UNION
    (SELECT DISTINCT
        p_pair.cid_1 cid_1, 
        p_pair.cid_2 cid_2, 
        p_pair.ensp_2 ensp_1, 
        p_pair.ensp_1 ensp_2
    FROM p_pair
    WHERE ensp_1 > ensp_2) );



