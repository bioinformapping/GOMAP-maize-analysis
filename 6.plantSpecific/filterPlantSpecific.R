source("R/gafTools.R")
source("R/oboTools.R")
source("R/getNonredDataset.R")
library("data.table")
library("argparse")


# getAllGo <- function(go_obo,x){
#   allGo = go_obo$ancestors[[x]]
#   if(!is.null(allGo)){
#     allGo = allGo
#   }else{
#     allGo = x
#   }
#   allGo
# }

parser <- ArgumentParser()

parser$add_argument("-i", "--input", type="character", default=TRUE,
  help="Input GOMAP gaf file", required=T)
parser$add_argument("-g", "--go_obo", type="character", default=TRUE,
  help="The go.obo file downloaded from geneontology.org", required=T)
parser$add_argument("-f", "--specific", type="character", default=TRUE,
  help="List of plant specific terms from the GOMAP-maize-analysis repo", required=T)
parser$add_argument("-o", "--output", type="character", default=TRUE,
  help="Output file to write plant-specific terms to", required=T)

args <- parser$parse_args()

go_obo <- check_obo_data(args$go_obo)
sppSpecTerms <- readLines(args$specific)
plantSpecificGO <- c(sppSpecTerms,c("GO:0005575","GO:0008150","GO:0003674"))

go_obo <- check_obo_data("data/go/go.obo")

getNRterms <- function(go_obo,data){
  out <- data[term_accession %in% minimal_set(go_obo,term_accession)]
  out
}

options(mc.cores = 8)

cat("Processing input file: ",args$input,'\n')
#input_gaf_data = read_gaf("data/clean/maize.PH207.UIUC_UMN-1.0.community.gaf")
input_gaf_data <- read_gaf(args$input)
input_gaf_data[,qualifier:=0]
input_gaf_data[,with:=0]
input_gaf_data[,db_object_name:=""]
input_gaf_data[,db_object_synonym:=""]
unique_terms = unique(input_gaf_data[,.(aspect,term_accession)])
all_terms = unique_terms[,list(all_term_accession=unique(go_obo$ancestors[[term_accession]])),by=list(aspect,term_accession)]
expand_gaf_data = unique(merge.data.table(input_gaf_data,all_terms,allow.cartesian = T,all.x = T))
expand_gaf_data[is.na(all_term_accession),all_term_accession:=term_accession]
expand_gaf_data = unique(expand_gaf_data[,-c("term_accession"),with=F])
colnames(expand_gaf_data)[colnames(expand_gaf_data)=="all_term_accession"] = "term_accession"
setcolorder(expand_gaf_data,gaf_cols)
today=format(Sys.Date(),"%Y%m%d")
expand_gaf_data[,date:=as.numeric(today)]
tmp_dt = data.table(plantSpecificGO)
colnames(tmp_dt) = "term_accession"
spp_gaf_data = unique(merge.data.table(tmp_dt,expand_gaf_data))
spp_list <- split(spp_gaf_data[,.(db_object_symbol,aspect,term_accession)],by=c("db_object_symbol","aspect"),keep.by = T)
spp_minimal_list <- mclapply(spp_list,getNRterms,go_obo=go_obo)
spp_minimal_dt <- unique(rbindlist(spp_minimal_list))
output_gaf_data <- unique(merge.data.table(spp_minimal_dt,spp_gaf_data,by=colnames(spp_minimal_dt),no.dups = T))
today=format(Sys.Date(),"%Y%m%d")
output_gaf_data[,date:=as.numeric(today)]
setcolorder(output_gaf_data,gaf_cols)
write_gaf(output_gaf_data,args$output)

