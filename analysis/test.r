library("data.table")
source("code/cafa_eval.r")
source("code/aigo_eval.r")

speciesSpecific <- fread("go/speciesSpecificGOTerms.txt")

speciesSpecific[`NCBITaxon:33090` != `NCBITaxon:3702`]
speciesSpecific[GOterm=="GO:0010229"]

all_go_file="go/go.obo"

#check and read/load the go.obo file
all_obo = check_obo_data(go_file)

plant_go_file

gaf_data <- read_gaf("annotations/goa_arabidopsis.gaf")


go_obo = check_obo_data("go/go.obo")

diff_terms <- setdiff(gaf_data[evidence_code!="IEA"]$term_accession,plant_obo$id)

diff_terms
go_obo$name[diff_terms]

setdiff(gaf_data[evidence_code!="IEA"]$term_accession,go_obo$id)
