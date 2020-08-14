----Training sets
CREATE TABLE chemical_interactions_v2.predictor_1 (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int,
);

CREATE TABLE chemical_interactions_v2.predictor_2 (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int
);

CREATE TABLE chemical_interactions_v2.predictor_1_high_scored (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int
);

CREATE TABLE chemical_interactions_v2.predictor_2_high_scored (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int
);

CREATE TABLE chemical_interactions_v2.predictor_1_independent_validation (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int
);

CREATE TABLE chemical_interactions_v2.predictor_2_independent_validation (
is_synergistic varchar,
ensp_1 int,
ensp_2 int,
cid_1 int,
cid_2 int
);



----Drugs
CREATE TABLE chemical_interactions_v2.all_tested_drugs (
stitch_pubchem_cid int, 
pubchem_cid_sub int,
approved varchar, 
pos_drugs integer,
atc_bc int,
side_effect_bc int,
shared_targets_with_pos_drug_bc int,
smile text,
atc varchar, 
side_effect varchar);


CREATE TABLE chemical_interactions_v2.drug_atc (
is_synergistic varchar,
cid_1 int,
cid_2 int,
atc_similarity real
);


CREATE TABLE chemical_interactions_v2.drug_side_effect (
is_synergistic varchar,
cid_1 int,
cid_2 int,
side_effect_similarity real
);



----Targets
CREATE TABLE chemical_interactions_v2.target_alias(
ensp int,
ensg int,
uniprot_id varchar,
entrez_kegg_id int,
gene_name varchar,
biogrid_id int
);



----FxnAss
CREATE TABLE chemical_interactions_v2.string_topo(
type varchar,
ensp_1 int,
ensp_2 int,
degree double precision,
clustering_coef double precision,
betweenness double precision,
closeness double precision,
neighbor_sharing double precision,
shortest_path_length double precision
);

CREATE TABLE chemical_interactions_v2.string_topo_eb(
type varchar,
ensp_1 int,
ensp_2 int,
eb_max double precision,
eb_min double precision,
eb_mean double precision,
eb_fraction double precision
);

CREATE TABLE chemical_interactions_v2.string_topo_nsp(
type varchar,
ensp_1 int,
ensp_2 int,
number_of_shortest_paths int);



----GeneEss
CREATE TABLE chemical_interactions_v2.gene_essentiality(
ensp int,
gene_essentiality smallint);



----GeneExp
CREATE TABLE chemical_interactions_v2.gene_coexpression(
ensp_1 int,
ensp_2 int,
coexpression double precision);



----GetInt
CREATE TABLE chemical_interactions_v2.getint_topo(
type varchar,
ensp_1 int,
ensp_2 int,
degree double precision,
clustering_coef double precision,
betweenness double precision,
closeness double precision,
neighbor_sharing double precision,
shortest_path_length double precision
);

CREATE TABLE chemical_interactions_v2.getint_topo_eb(
type varchar,
ensp_1 int,
ensp_2 int,
eb_max double precision,
eb_min double precision,
eb_mean double precision,
eb_fraction double precision
);

CREATE TABLE chemical_interactions_v2.getint_topo_nsp(
type varchar,
ensp_1 int,
ensp_2 int,
number_of_shortest_paths int);



----PPI
CREATE TABLE chemical_interactions_v2.biogrid_topo(
type varchar,
ensp_1 int,
ensp_2 int,
degree double precision,
clustering_coef double precision,
betweenness double precision,
closeness double precision,
neighbor_sharing double precision,
shortest_path_length double precision
);

CREATE TABLE chemical_interactions_v2.biogrid_topo_eb(
type varchar,
ensp_1 int,
ensp_2 int,
eb_max double precision,
eb_min double precision,
eb_mean double precision,
eb_fraction double precision
);

CREATE TABLE chemical_interactions_v2.biogrid_topo_nsp(
type varchar,
ensp_1 int,
ensp_2 int,
number_of_shortest_paths int);



----GO
CREATE TABLE chemical_interactions_v2.go_all(
ensp_1 int,
ensp_2 int,
go_all_sem_sim double precision);

CREATE TABLE chemical_interactions_v2.go_bp(
ensp_1 int,
ensp_2 int,
go_bp_sem_sim double precision);

CREATE TABLE chemical_interactions_v2.go_cc(
ensp_1 int,
ensp_2 int,
go_cc_sem_sim double precision);

CREATE TABLE chemical_interactions_v2.go_mf(
ensp_1 int,
ensp_2 int,
go_mf_sem_sim double precision);

CREATE TABLE chemical_interactions_v2.go_slim(
ensp_1 int,
ensp_2 int,
go_slim_sem_sim double precision);



---- Phylo
CREATE TABLE chemical_interactions_v2.phylo(
ensp_1 int,
ensp_2 int,
phylogenic_similarity double precision);




