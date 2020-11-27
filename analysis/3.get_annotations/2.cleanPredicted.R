source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
library("data.table")
library("yaml")

config = read_yaml("config.yml")
go_obo <- check_obo_data("data/go/go.obo")

predData_list <- lapply(config$predicted,function(x){
  inputGaf = file.path("data","predicted",x)
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,qualifier:=0]
  input_gaf_data[,with:=0]
  input_gaf_data[,db_object_name:=""]
  input_gaf_data[,db_object_synonym:=""]
  gaf_cols[gaf_cols!="term_accession"]
  output_gaf_data <- input_gaf_data
  today=format(Sys.Date(),"%Y%m%d")
  output_gaf_data[,date:=as.numeric(today)]
  output_gaf_data[term_accession %in% names(go_obo$alt_conv),
                  term_accession:=go_obo$alt_conv[term_accession]]
  outputGaf = file.path("data","clean",x)
  write_gaf(output_gaf_data,outputGaf)
  NULL
})
