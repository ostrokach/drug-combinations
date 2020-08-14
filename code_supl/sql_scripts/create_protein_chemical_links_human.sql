DROP TABLE stitch.protein_chemical_links_human;

CREATE TABLE stitch.protein_chemical_links_human (
chemical varchar,
cid_full int,
cid int,
protein varchar,
ensp int,
combined_score int,
rank int);

INSERT INTO stitch.protein_chemical_links_human (chemical, protein, combined_score)
SELECT 
chemical, 
protein, 
combined_score
FROM stitch.protein_chemical_links
WHERE protein LIKE '%ENSP0%';

UPDATE stitch.protein_chemical_links_human
SET cid_full = substr(chemical, 4)::int;

UPDATE stitch.protein_chemical_links_human
SET cid = substr(chemical, 5)::int;

UPDATE stitch.protein_chemical_links_human
SET ensp = regexp_replace(protein, '.+\.ENSP', '')::int;

CREATE INDEX protein_chemical_links_human_cid_ensp_idx
  ON stitch.protein_chemical_links_human
  USING btree
  (cid , ensp );

CREATE INDEX protein_chemical_links_human_cid_idx
  ON stitch.protein_chemical_links_human
  USING btree
  (cid );
  
CREATE INDEX protein_chemical_links_human_ensp_idx
  ON stitch.protein_chemical_links_human
  USING btree
  (ensp );


DROP TABLE stitch.protein_chemical_links_human_nostereo;

CREATE TABLE stitch.protein_chemical_links_human_nostereo (
cid int,
ensp int,
combined_score int,
rank int);

INSERT INTO stitch.protein_chemical_links_human_nostereo (cid, ensp, combined_score)
SELECT
cid, 
ensp, 
MAX(combined_score)
FROM stitch.protein_chemical_links_human
GROUP BY cid, ensp;

CREATE UNIQUE INDEX ON stitch.protein_chemical_links_human_nostereo (cid, ensp);

CREATE INDEX ON stitch.protein_chemical_links_human_nostereo (cid);

CREATE INDEX ON stitch.protein_chemical_links_human_nostereo (ensp);

WITH pcl_temp AS (
  SELECT 
    cid,
    ensp,
    RANK() OVER (PARTITION BY cid ORDER BY combined_score DESC) rk
  FROM stitch.protein_chemical_links_human_nostereo)
UPDATE stitch.protein_chemical_links_human_nostereo pcl
SET rank = pcl_temp.rk
FROM pcl_temp
WHERE pcl.cid = pcl_temp.cid
AND pcl.ensp = pcl_temp.ensp;

CREATE INDEX ON stitch.protein_chemical_links_human_nostereo (rank);
