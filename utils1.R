library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(KEGGREST)
library(Rdisop)
library(ontologyIndex)
source("functions.R")

utils1.blast_all2all = function()
{
  system("diamond makedb --in data/db/UNIPROT/ecr_proteins.faa -d data/db/UNIPROT/ecr_proteins")
  system("diamond blastp --sensitive --max-target-seqs 10000 -d data/db/UNIPROT/ecr_proteins -q data/db/UNIPROT/csh_proteins.faa -o data/db/UNIPROT/csh_ecr.m8 -f 6 qseqid sseqid qlen slen score bitscore length")
}

utils1.collect_kegg = function()
{
  uniprot2kegg = readr::read_tsv("data/db/uniprot2kegg.tsv", col_names=T, na="")
  kegg.species2org = readr::read_delim("data/db/KEGG/kegg_species2org.tsv", "\t")
  kegg.gene2enzyme = do.call(rbind, lapply(kegg.species2org$KEGG_ORG, function(org) {
    response = KEGGREST::keggLink(org, target="enzyme")
    data.frame(KEGG_ORG=org, KEGG_GENE=names(response), KEGG_EC=response)
  })) %>% dplyr::left_join(uniprot2kegg, by=c("KEGG_GENE"="kegg_protein_id")) %>%
    dplyr::rename(KEGG_UNIPROT="uniprot_id")
  readr::write_delim(kegg.gene2enzyme, "data/db/KEGG/kegg.gene2enzyme.tsv", delim="\t", na="")
  
  kegg.enzyme2cpd = stack(KEGGREST::keggLink("enzyme", "cpd")) %>% dplyr::select(KEGG_CPD=ind, KEGG_EC=values)
  
  kegg.compound_ids = unique(as.character(kegg.enzyme2cpd$KEGG_CPD))
  kegg.compounds = data.frame()
  for(ids in chunk(kegg.compound_ids, 10)) {
    response = KEGGREST::keggGet(dbentries=ids)
    if(length(response) < length(ids)) {
      asd()
    }
    for(i in 1:length(response)) {
      kegg.compounds.i = with(response[[i]], data.frame(
        KEGG_CPD=ENTRY, 
        KEGG_FORMULA=ifelse("FORMULA" %in% names(response[[i]]), FORMULA, NA_character_), 
        KEGG_EXACT_MASS=ifelse("EXACT_MASS" %in% names(response[[i]]), EXACT_MASS, NA_real_)))
      kegg.compounds = rbind(kegg.compounds, kegg.compounds.i)
    }
  }
  readr::write_delim(kegg.compounds, "data/db/KEGG/kegg_compounds.tsv", delim="\t", na="")
  
  kegg.enzyme2cpd = kegg.enzyme2cpd %>%
    dplyr::inner_join(kegg.compounds %>% dplyr::mutate(KEGG_CPD=paste0("cpd:", KEGG_CPD)), by="KEGG_CPD")
  readr::write_delim(kegg.enzyme2cpd, "data/db/KEGG/kegg_enzyme2cpd.tsv", delim="\t", na="")
  
  kegg.species2cpd = kegg.enzyme2cpd %>%
    dplyr::inner_join(kegg.gene2enzyme, by="KEGG_EC")%>%
    dplyr::inner_join(kegg.species2org, by="KEGG_ORG")  %>%
    dplyr::distinct(Species, Species2, KEGG_ORG, KEGG_GENE, KEGG_UNIPROT, KEGG_CPD, KEGG_FORMULA, KEGG_EXACT_MASS)
  readr::write_delim(kegg.species2cpd, "data/db/KEGG/kegg_species2cpd.tsv", delim="\t", na="")
  
  
  kegg.species2pathways = do.call(rbind, lapply(kegg.species2org$KEGG_ORG, function(org) {
    stack(KEGGREST::keggLink("pathway", org)) %>%
      dplyr::mutate(KEGG_ORG=org, values=gsub(org, "map", values)) %>%
      dplyr::select(KEGG_ORG, KEGG_PATHWAY=values) %>%
      dplyr::distinct()
  }))
  
  kegg.pathways2cpd = stack(KEGGREST::keggLink("pathway", "compound")) %>%
    dplyr::rename(KEGG_PATHWAY="values", KEGG_CPD="ind") %>%
    dplyr::inner_join(kegg.species2cpd %>% dplyr::distinct(KEGG_ORG, KEGG_CPD), by="KEGG_CPD") %>%
    dplyr::inner_join(kegg.compounds %>% dplyr::mutate(KEGG_CPD=paste0("cpd:", KEGG_CPD)), by="KEGG_CPD") %>%
    dplyr::inner_join(kegg.species2pathways, by=c("KEGG_ORG", "KEGG_PATHWAY"))
  readr::write_delim(kegg.pathways2cpd, "data/db/KEGG/kegg_pathway2cpd.tsv", delim="\t", na="")
  
  kegg.pathways = stack(KEGGREST::keggList("pathway")) %>%
    dplyr::rename(KEGG_PATHWAY_NAME="values", KEGG_PATHWAY="ind")
  readr::write_delim(kegg.pathways, "data/db/KEGG/kegg_pathways.tsv", delim="\t", na="")
  
  kegg.pathway2gene = stack(KEGGREST::keggLink("pathway", "enzyme")) %>%
    dplyr::select(KEGG_EC=ind, KEGG_PATHWAY=values) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(kegg.gene2enzyme, by="KEGG_EC") %>%
    dplyr::mutate(KEGG_PATHWAY=gsub(":ec", ":map", KEGG_PATHWAY)) %>%
    dplyr::select(KEGG_ORG, KEGG_UNIPROT, KEGG_GENE, KEGG_PATHWAY, KEGG_EC) %>%
    dplyr::distinct()
  readr::write_delim(kegg.pathway2gene, "data/db/KEGG/kegg.pathway2gene.tsv", delim="\t", na="")
}

utils1.collect_pathways2peaks = function()
{
  ppm_extended = 100e-6
  ppm = 10e-6
  electron_mass = 0.00054858026 
  adducts = readr::read_delim("data/db/adducts.tsv", "\t", comment="#") %>%
    tidyr::separate(adduct_formula, into=c("adduct_formula_add", "adduct_formula_sub"), sep="-") %>% 
    dplyr::filter(adduct_name %in% c("[M-2H]2-", "[M-H]-"))
  kegg.pathways2formula = readr::read_delim("data/db/KEGG/kegg_pathway2cpd.tsv", delim="\t") %>%
    dplyr::filter(grepl("map00230", KEGG_PATHWAY) & !is.na(KEGG_FORMULA) & !grepl("R", KEGG_FORMULA) & KEGG_EXACT_MASS > 50) %>%
    dplyr::distinct(KEGG_PATHWAY, KEGG_FORMULA, KEGG_CPD)
  kegg.pathways2peaks = kegg.pathways2formula %>%
    tidyr::crossing(adducts) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z) {
      z.molecule = Rdisop::getMolecule(z$KEGG_FORMULA)
      if(z$adduct_nmol > 1) {
        for(i in 2:z$adduct_nmol) { z.molecule = Rdisop::addMolecules(z.molecule$formula, z$formula) }
      }
      
      if(!is.na(z$adduct_formula_add) && z$adduct_formula_add!="") {
        z.molecule = Rdisop::addMolecules(z.molecule$formula, z$adduct_formula_add)
      }
      if(!is.na(z$adduct_formula_sub) && z$adduct_formula_sub!="") {
        z.molecule = Rdisop::subMolecules(z.molecule$formula, z$adduct_formula_sub)
      }
      ion_mz = (z.molecule$isotopes[[1]][1,] - z$adduct_charge*electron_mass)/abs(z$adduct_charge)
      cbind(data.frame(z), abundance=z.molecule$isotopes[[1]][2,], mz=ion_mz)
    })(.)) %>%
    dplyr::mutate(lb=mz*(1-ppm), ub=mz*(1+ppm), lb_extended=mz*(1-ppm_extended), ub_extended=mz*(1+ppm_extended))
  
  readr::write_delim(kegg.pathways2peaks, "data/db/KEGG/kegg_pathway2peaks.tsv", delim="\t", na="", col_names=T)
}

utils1.collect_tigrfam_info = function()
{
  tigrfam.info = data.frame()
  for(file in list.files("data/db/TIGRFAM/info", pattern="INFO$", full.names=T)) {
    f = read_lines(file)
    f = stringr::str_match(f, "([^ ]+) +(.*)")
    vars = f[,3]
    names(vars) = f[,2]
    
    cutoff.trusted = as.numeric(strsplit(vars["TC"], " ")[[1]][1])
    cutoff.noise = as.numeric(strsplit(vars["NC"], " ")[[1]][1])
    ec = vars["EC"]
    tigrfam.info.f = data.frame(annotation=vars["AC"], cutoff.trusted=cutoff.trusted, cutoff.noise=cutoff.noise, type=vars["IT"], ec=ec)
    tigrfam.info = rbind(tigrfam.info, tigrfam.info.f)
  }
  
  readr::write_delim(tigrfam.info, path="data/db/TIGRFAM/info.tsv", delim="\t", na="", col_names=T)
}

utils1.collect_goterms = function()
{
  go_obo = ontologyIndex::get_OBO("http://purl.obolibrary.org/obo/go/snapshot/go.obo", propagate_relationships = "is_a", extract_tags = "minimal")
  go_ontology = sapply(go_obo$ancestors, function(z) dplyr::case_when(any(z=="GO:0008150") ~ "biological_process", any(z=="GO:0003674") ~ "molecular_function", any(z=="GO:0005575") ~ "cellular_component", T ~ "obsolete"))
  go_terms = data.frame(go_term=go_obo$id, go_description=go_obo$name, go_ontology=go_ontology) %>%
    dplyr::filter(grepl("^GO:", go_term))
  
  readr::write_delim(go_terms, path="data/db/GO/go_terms.tsv", delim="\t", na="", col_names=T)
}

utils1.collect_families = function() 
{
  uniprot2family = data.frame()
  for(file in list.files("data/db/UNIPROT", pattern="proteinids2refseq.list$", full.names=T)) {
    org = gsub("_.*", "", basename(file))
    protein2refseq = readr::read_delim(file, "\t") %>% dplyr::select(protein_id="From", accession="To")
    refseq2eggnog = readr::read_delim(paste0("data/db/EGGNOG/", gsub("_proteinids2refseq.list", "_eggnog.tsv", basename(file))), "\t", locale=locale(asciify=T), comment="#", col_names=c("accession", "seed_eggnog_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "best_tax_level", "name", "go", "ec", "kegg_ko", "kegg_pathway", "kegg_module", "kegg_reaction", "kegg_rclass", "brite", "kegg_tc", "CAZy", "bigg_reaction", "taxonomy", "eggnog", "smthng", "cog", "description"))
    refseq2tigrfam = readr::read_delim(paste0("data/db/TIGRFAM/", gsub("_proteinids2refseq.list", "_tigrfam.tsv", basename(file))), "\t", locale=locale(asciify=T), comment="#", col_names=c("tigrfam", "target", "accession", "accession2", "evalue", "score", "bias", "evalue_bestdomain", "score_bestdomain", "bias_bestdomain", "domains_expected", "regions_number", "multidomain_number", "overlaps_number", "envelopes_number", "domains_number", "domains_reported_number", "domains_significant_number", "description"))
  
    uniprot2family.eggnog = refseq2eggnog %>%
      dplyr::select(accession, score=seed_ortholog_score, evalue=seed_ortholog_evalue, cluster=eggnog, kegg_ko, kegg_pathway, go, description) %>%
      reshape2::melt(id.vars=c("accession", "score", "evalue", "description"), variable.name="annotation_source", value.name="annotation") %>%
      dplyr::filter(!is.na(annotation)) %>%
      tidyr::separate_rows(annotation, sep=",") %>%
      dplyr::filter(!is.na(annotation)) %>%
      dplyr::mutate(annotation_source=paste0("eggnog_v2.", annotation_source)) %>%
      data.frame() %>%
      dplyr::distinct(accession, annotation_source, annotation,  .keep_all=T)
    uniprot2family.tigrfam = refseq2tigrfam %>%
      dplyr::select(accession, score, evalue, family=tigrfam, description) %>%
      reshape2::melt(id.vars=c("accession", "score", "evalue", "description"), variable.name="annotation_source", value.name="annotation") %>%
      dplyr::filter(!is.na(annotation)) %>%
      tidyr::separate_rows(annotation, sep=",") %>%
      dplyr::filter(!is.na(annotation)) %>%
      dplyr::mutate(annotation_source=paste0("tigrfam.", annotation_source)) %>%
      data.frame() %>%
      dplyr::distinct(accession, annotation_source, annotation,  .keep_all=T)
    uniprot2family.f = rbind(uniprot2family.eggnog, uniprot2family.tigrfam) %>% 
      dplyr::mutate(org=org, annotation_id=paste0(annotation_source, "/", annotation)) %>%
      dplyr::inner_join(protein2refseq, by="accession")
    uniprot2family = rbind(uniprot2family, uniprot2family.f)
  }
  
  write.table(uniprot2family, "data/db/uniprot2family.tsv", na="", sep="\t", quote=F, row.names=F, col.names=T)
}