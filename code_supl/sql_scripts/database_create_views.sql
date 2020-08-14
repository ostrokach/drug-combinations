-- Predictor 0

COPY (

SELECT

predictor.*,
drug_atc_similarity.atc_similarity atc_similarity,
drug_chemical_similarity.morganfingerprintr2_tanimoto chemical_similarity,
drug_side_effect_similarity.side_effect_similarity side_effect_similarity,

FROM

chemical_interactions_v2.predictor_1_all_unused_pairs_3 predictor
JOIN chemical_interactions_v2.drug_atc_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_chemical_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_side_effect_similarity USING (cid_1, cid_2)

) TO '/tmp/predictor_0_all_unused_pairs_3.tsv' WITH CSV DELIMITER E'\t' HEADER;


-- Predictor 1

COPY (

SELECT

predictor.*,
drug_atc_similarity.atc_similarity atc_similarity,
drug_chemical_similarity.morganfingerprintr2_tanimoto chemical_similarity,
drug_side_effect_similarity.side_effect_similarity side_effect_similarity,
biogrid_topo.shortest_path_length biogrid_shortest_path_length,
biogrid_topo_eb.eb_max biogrid_eb_max,
gene_coexpression.coexpression gene_coexpression,
gene_ess_1.gene_essentiality gene_essentiality_1,
gene_ess_2.gene_essentiality gene_essentiality_2,
getint_topo.shortest_path_length getint_shortest_path_length,
getint_topo_eb.eb_max getint_eb_max,
go_all.go_all_sem_sim go_all_sem_sim,
go_bp.go_bp_sem_sim go_bp_sem_sim,
go_cc.go_cc_sem_sim go_cc_sem_sim,
go_mf.go_mf_sem_sim go_mf_sem_sim,
phylo.phylogenic_similarity phylogenic_similarity,
string_topo.shortest_path_length string_shortest_path_length,
string_topo_eb.eb_max string_eb_max

FROM

chemical_interactions_v2.predictor_1_all_unused_pairs_3 predictor
JOIN chemical_interactions_v2.drug_atc_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_chemical_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.drug_side_effect_similarity USING (cid_1, cid_2)
JOIN chemical_interactions_v2.biogrid_topo USING (ensp_1, ensp_2)
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

) TO '/tmp/predictor_1_all_unused_pairs_3.tsv' WITH CSV DELIMITER E'\t' HEADER;



-- Predictor 2

COPY (

SELECT

predictor.*,
biogrid_topo.shortest_path_length biogrid_shortest_path_length,
biogrid_topo_eb.eb_max biogrid_eb_max,
gene_coexpression.coexpression gene_coexpression,
gene_ess_1.gene_essentiality gene_essentiality_1,
gene_ess_2.gene_essentiality gene_essentiality_2,
getint_topo.shortest_path_length getint_shortest_path_length,
getint_topo_eb.eb_max getint_eb_max,
go_all.go_all_sem_sim go_all_sem_sim,
go_bp.go_bp_sem_sim go_bp_sem_sim,
go_cc.go_cc_sem_sim go_cc_sem_sim,
go_mf.go_mf_sem_sim go_mf_sem_sim,
phylo.phylogenic_similarity phylogenic_similarity,
string_topo.shortest_path_length string_shortest_path_length,
string_topo_eb.eb_max string_eb_max

FROM

chemical_interactions_v2.predictor_2_all_unused_pairs_3
JOIN chemical_interactions_v2.biogrid_topo USING (ensp_1, ensp_2)
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

) TO '/tmp/predictor_2_all_unused_pairs_3.tsv' WITH CSV DELIMITER E'\t' HEADER;

