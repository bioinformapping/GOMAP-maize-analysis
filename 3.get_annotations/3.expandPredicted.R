source("R/gafTools.R")
source("R/oboTools.R")
library("data.table")
library("yaml")

config = read_yaml("config.yaml")
go_obo <- check_obo_data("data/go/go.obo")

getAllGo <- function(go_obo,x){
  allGo = go_obo$ancestors[[x]]
  if(!is.null(allGo)){
    allGo = allGo
  }else{
    allGo = x
  }
  allGo
}


predData_list <- lapply(config$predicted,function(x){
  print(x["inbred"])
  inputGaf = file.path("data","clean",x$file)
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,qualifier:=0]
  input_gaf_data[,with:=0]
  input_gaf_data[,db_object_name:=""]
  input_gaf_data[,db_object_synonym:=""]
  unique_terms = unique(input_gaf_data[,.(aspect,term_accession)])
  all_terms = unique_terms[,list(all_term_accession=unique(go_obo$ancestors[[term_accession]])),by=list(aspect,term_accession)]
  expand_gaf_data = unique(merge.data.table(input_gaf_data,all_terms,allow.cartesian = T,all.x = T))
  expand_gaf_data[is.na(all_term_accession),all_term_accession:=term_accession]
  output_gaf_data = unique(expand_gaf_data[,-c("term_accession"),with=F])
  colnames(output_gaf_data)[colnames(output_gaf_data)=="all_term_accession"] = "term_accession"
  setcolorder(output_gaf_data,gaf_cols)
  today=format(Sys.Date(),"%Y%m%d")
  output_gaf_data[,date:=as.numeric(today)]
  outputGaf = file.path("data","expanded",x$file)
  write_gaf(output_gaf_data,outputGaf)
  NULL
})

