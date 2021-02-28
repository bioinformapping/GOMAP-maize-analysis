source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
library("data.table")

sppSpecTerms <- fread("data/go/speciesSpecificGOTerms.txt")
plantSpecificGO <- sppSpecTerms[`NCBITaxon:33090`==1]$GOterm
plantSpecificGO <- c(plantSpecificGO,c("GO:0005575","GO:0008150","GO:0003674"))

inbredLines <- c("B73v3","B73v4","Mo17","PH207","W22")
allGoFiles <- dir("data/7.expanded",full.names = T)
names(allGoFiles) = inbredLines

go_obo <- check_obo_data("data/go/go.obo")mess

getNRterms <- function(terms,sppSpecific){
  sppTerms = minimal_set(go_obo,terms[which(sppSpecific==T)])
  out <- list(term_accession=sppTerms)
  out
}

allGoData_list <- lapply(inbredLines,function(x){
  inputGaf <- allGoFiles[[x]]
  cat("Processing input file",inputGaf,'\n')
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,sppSpecific:=ifelse(term_accession %in% plantSpecificGO,TRUE,FALSE)]
  tmp_gaf_data <- input_gaf_data[,getNRterms(term_accession,sppSpecific),by=list(db_object_symbol,aspect)]
  #print(tmp_gaf_data)
  output_gaf_data <- merge.data.table(tmp_gaf_data,unique(input_gaf_data[,-c("sppSpecific","term_accession"),with=F]),no.dups = T)
  #print(output_gaf_data)
  today=format(Sys.Date(),"%d%m%Y")
  output_gaf_data[,date:=as.numeric(today)]
  outputGaf = gsub(".gz$","",gsub("7.expanded","8.specific",inputGaf))
  cat("Writing output file",outputGaf,"\n")
  write_gaf(output_gaf_data,outputGaf)
  input_gaf_data
})

