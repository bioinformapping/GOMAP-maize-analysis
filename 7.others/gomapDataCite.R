library("data.table")
library("bib2df")


dataset_dt = fread("data/citations/gomap_datasets.tsv",sep = "\t",header = T)
dataset_bib = data.table(bib2df("data/citations/gomap_datasets.bib"))

dataset_dt[,c("first","last"):=tstrsplit(Person," ")]

dataset_dt[,author:=paste0(last,",",first)]
paste0(dataset_dt$last,"_",dataset_dt$Rel,"_",gsub(" +","_",dataset_dt$Line))

author_list = lapply(dataset_dt$author, function(x){
  c(x,"Larence-Dill, Carolyn")
})
dataset_bib$AUTHOR = author_list
dataset_bib$doi = dataset_dt$DOI
dataset_bib$URL = paste0("http://dx.doi.org/",dataset_dt$DOI)
#dataset_bib$YEAR = gsub("},","",dataset_bib$YEAR,fixed = T)

dataset_bib$ABSTRACT = gsub("_","\\_",dataset_bib$NOTE,fixed=T)
dataset_bib$TITLE = gsub("_","\\_",dataset_bib$TITLE,fixed=T)
dataset_bib$BIBTEXKEY = paste0(dataset_dt$last,"_",dataset_dt$Rel,"_",gsub(" +","_",dataset_dt$Line))

dataset_bib[1]
bib2df::df2bib(dataset_bib[,-c("NOTE"),with=F],"data/citations/test.bib")


test_bib = bib2df("data/citations/gomap_datasets.bib")

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
