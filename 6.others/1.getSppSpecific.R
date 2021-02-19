source("R/gafTools.R")
source("R/oboTools.R")
source("R/getNonredDataset.R")
library("data.table")
library("yaml")

sppSpecTerms <- fread("data/go/speciesSpecificGOTerms.txt")
plantSpecificGO <- sppSpecTerms[`NCBITaxon:33090`==1]$GOterm
plantSpecificGO <- c(plantSpecificGO,c("GO:0005575","GO:0008150","GO:0003674"))

config = read_yaml("config.yaml")
go_obo <- check_obo_data("data/go/go.obo")

getNRterms <- function(go_obo,data){
  out <- data[term_accession %in% minimal_set(go_obo,term_accession)]
  out
}

options(mc.cores = 8)
#x=config$predicted[[1]]
allGoData_list <- lapply(config$predicted,function(x){
  inputGaf <- file.path("data","expanded",x["file"])
  cat("Processing input file: ",inputGaf,'\n')
  input_gaf_data <- read_gaf(inputGaf)
  tmp_dt = data.table(plantSpecificGO)
  colnames(tmp_dt) = "term_accession"
  spp_gaf_data = merge.data.table(tmp_dt,input_gaf_data)
  spp_list <- split(spp_gaf_data[,.(db_object_symbol,aspect,term_accession)],by=c("db_object_symbol","aspect"),keep.by = T)
  spp_mininal_list <- mclapply(spp_list,getNRterms,go_obo=go_obo)  
  spp_mininal_dt <- rbindlist(spp_mininal_list)
  output_gaf_data <- merge.data.table(spp_gaf_data,spp_mininal_dt,no.dups = T)
  today=format(Sys.Date(),"%Y%m%d")
  output_gaf_data[,date:=as.numeric(today)]
  outputGaf = file.path("data","plant-specific",x["file"])
  setcolorder(output_gaf_data,gaf_cols)
  write_gaf(output_gaf_data,outputGaf)
  input_gaf_data
})
