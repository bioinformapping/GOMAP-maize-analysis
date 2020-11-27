source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
library("data.table")
library("yaml")

config = read_yaml("config.yml")
go_obo <- check_obo_data("data/go/go.obo")

getAllGo <- function(go_obo,x){
  allGo = unique(unlist(go_obo$ancestors[x]))
  if(!is.null(allGo)){
    allGo = allGo
  }else{
    allGo = x
  }
  allGo
}

predData_list <- mclapply(config$predicted,function(x){
  print(x)
  inputGaf = file.path("data","clean",x)
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,qualifier:=0]
  input_gaf_data[,with:=0]
  input_gaf_data[,db_object_name:=""]
  input_gaf_data[,db_object_synonym:=""]
  output_gaf_data <- unique(input_gaf_data[,list(term_accession=getAllGo(go_obo,term_accession)),by=eval(gaf_cols[!gaf_cols %in% "term_accession"])])
  today=format(Sys.Date(),"%Y%m%d")
  output_gaf_data[,date:=as.numeric(today)]
  outputGaf = file.path("data","expanded",x)
  write_gaf(output_gaf_data,outputGaf)
  outputGaf
},mc.cores = 4)

