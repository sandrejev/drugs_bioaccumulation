library(ggplot2)

#
# t-test
#
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) {
    # welch-satterthwaite equation
    se = sqrt( (s1^2/n1) + (s2^2/n2) )
    df = ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else {
    # pooled standard deviation, scaled by the sample sizes
    se = sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df = n1+n2-2
  }      
  
  t = (m1-m2-m0)/se 
  dat = data.frame("Difference of means"=m1-m2, "Std Error"=se, "t"=t, "p_value"=2*pt(-abs(t),df))
  return(dat) 
}

chembl.molecules = function(ids) {
  ids_str = paste0(na.omit(ids), collapse=";")
  r = httr::GET(paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/set/", ids_str, "?format=xml"))
  xml = xml2::read_xml(httr::content(r, "text"))
  
  nodes = xml2::xml_find_all(xml, "//*[local-name()='molecule']")
  chembl_id = sapply(nodes, function(node) { v = xml2::xml_text(xml2::xml_find_all(node, "./molecule_chembl_id/text()")); if(length(v)==0) { return(NA) } else { return(v)} }) 
  chembl_name = sapply(nodes, function(node) { name = xml2::xml_text(xml2::xml_find_all(node, "./molecule_synonyms/synonym[1]/molecule_synonym/text()")); if(length(name)==0) { return(NA) } else { return(name)} })
  chembl_smiles = sapply(nodes, function(node) { v = xml2::xml_text(xml2::xml_find_all(node, "./molecule_structures/canonical_smiles/text()")); if(length(v)==0) { return(NA) } else { return(v)} })
  chembl_inchi = sapply(nodes, function(node) { v = xml2::xml_text(xml2::xml_find_all(node, "./molecule_structures/standard_inchi/text()")); if(length(v)==0) { return(NA) } else { return(v)} })
  
  
  data.frame(chembl_id=chembl_id, chembl_name=chembl_name, chembl_smiles=chembl_smiles, chembl_inchi=chembl_inchi, stringsAsFactors=F) %>%
    dplyr::filter(!is.na(chembl_id))
}

metlin.get = function(id) {
  for(r in list.files("/home/andrejev/Downloads/metlin", full.names = T)) {
    html = xml2::read_html(r)
    xml2::xml_find_all( html, "//*[contains(@href, 'pubchem.ncbi.nlm.nih.gov/compound')]")
    xml2::write_xml(html, file="test.xml")
    xml2::xml_find_all( html, "//a")
  }
}

metaboanalyst.mapper = function(ids, type="metlin") {
  x[1:2,]
  ids_str = paste0(ids, collapse=";")
  r = httr::POST(paste0("http://api.xialab.ca/mapcompounds"), httr::add_headers(`Content-Type`="application/json", `cache-control`="no-cache"), body=paste0('{"queryList": "',ids_str, '", "inputType": "', type, '"}'))
  response = httr::content(r, "text")
  json = jsonlite::parse_json(response)
  df = do.call(rbind, lapply(json, data.frame)) %>%
    dplyr::mutate_all(function(z) ifelse(z=="NA" | z=="-", NA_character_, as.character(z))) 
  df[!is.na(df[,paste0(type, "_id")]),]
}

kegg.mapping = function(dest) {
  r = httr::GET(paste0("http://rest.genome.jp/link/", dest, "/compound"))
  mappings = readr::read_tsv(httr::content(r, "text"), col_names = F) %>%
    dplyr::select(src=X1, dest=X2) %>%
    dplyr::mutate(src=gsub("^[^:]+:", "", src), dest=gsub("^[^:]+:", "", dest))
  
  mappings
}

unichem.mapping = function(src, dest) {
  url = paste0("https://www.ebi.ac.uk/unichem/rest/mapping/", dest, "/", src)
  print(url)
  r = httr::GET(url, httr::add_headers(`Accept`="text/xml", `Content-Type`="text/xml"))
  response = gsub(paste0(src, "="), "src=", gsub(paste0(dest, "="), "dest=", httr::content(r, "text")))
  xml = xml2::read_xml(response)
  mappings = do.call(rbind, lapply(xml2::xml_find_all(xml, "//*[local-name()='data']"), xml2::xml_attrs))
  as.data.frame(mappings, stringsAsFactors=F)
}

myTheme <- list(theme_minimal(base_size = 15), 
                theme(panel.background = element_rect(color = "grey90"), 
                      panel.grid.minor = element_line(colour = "white"), 
                      panel.grid.major = element_line(colour = "grey95")),
                geom_vline(xintercept = 0, color = "grey50") , 
                geom_hline(yintercept = 0, color = "grey50"),
                theme(axis.title.y = element_text(size = rel(0.8), angle = 90)),
                theme(axis.title.x = element_text(size = rel(0.8), angle = 00)),
                theme(legend.title = element_text(size = rel(0.8))),
                theme(legend.text = element_text(size = rel(0.8))))



log.na = function(x) log10(ifelse(x>0, x, NA))
log.zero = function(x) log10(ifelse(x>0, x, 1))
glog2 = function(x) ((asinh(x)-log(2))/log(2))


strset = function(x, f=T) {
  x = sort(unique(na.omit(x[f])))
  if(length(x)==0) return(NA_character_)
  paste0(x, collapse="|")
}

chunk = function(x,n) {
  split(x, rep(1:ceiling(length(x)/n), each=n)[1:length(x)])
}