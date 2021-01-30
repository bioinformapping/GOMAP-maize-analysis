suppressMessages({
  source("R/oboTools.R")
  source("R/gafTools.R")
  source("R/getNonredDataset.R")
  library("data.table")
})

go_obo <- check_obo_data("data/go/go.obo")
go2aspect <- list(biological_process="P",cellular_component="C",molecular_function="F")

# Processing the PH207 Annotations

print("Processing PH207")

ph207 <- fread("data/community/TPC2016-00353-LSBR2-A_Supplemental_Dataset_1.txt")
ph207Gaf <- ph207[,.(PH207_gene,GO_ID)]
ph207Gaf <- ph207[GO_ID!="N/A",list(term_accession=unlist(strsplit(GO_ID,";"))),by=list(PH207_gene)]
colnames(ph207Gaf)[1] <- "db_object_symbol"
ph207Gaf=ph207Gaf[,list(term_accession=goterm2alt_id(term_accession,go_obo)),by=db_object_symbol]
ph207Gaf[,db_object_id:=db_object_symbol]
ph207Gaf[,db:="MaizeGDB"]
ph207Gaf[,qualifier:=""]
ph207Gaf[,db_reference:=""]
ph207Gaf[,evidence_code:="IEA"]
ph207Gaf[,with:=""]
ph207Gaf[,db_object_name:=""]
ph207Gaf[,db_object_synonym:=""]
ph207Gaf[,db_object_type:="protein"]
ph207Gaf[,taxon:="taxon:4577"]
ph207Gaf[,date:="20161106"]
ph207Gaf[,assigned_by:="Hirsch, et al 2016"]
ph207Gaf[,c(setdiff(gaf_cols,colnames(ph207Gaf))):=""]
ph207Gaf[,aspect:=goterm2aspect(ph207Gaf$term_accession,go_obo)]
ph207Gaf_nonred <- rm_gaf_red(ph207Gaf,go_obo)
write_gaf(ph207Gaf_nonred,"data/predicted/maize.PH207.UIUC_UMN-1.0.community.gaf")

# Processing the Gramene B73v4 Annotations

print("Processing Gramene B73v4")

b73v4_gramene <- fread("data/community/zmays_refgeb_v4_gramene_go.txt")
b73v4_gramene_gaf <- unique(b73v4_gramene[grep("Zm00001d",`Gene stable ID`),c("Gene stable ID","GO term accession","GO term evidence code")])
colnames(b73v4_gramene_gaf) <- c("db_object_symbol","go_id","evidence_code")

b73v4_gramene_gaf=b73v4_gramene_gaf[,list(term_accession=goterm2alt_id(go_id,go_obo)),by=list(db_object_symbol,evidence_code)]
b73v4_gramene_gaf[,aspect:=goterm2aspect(term_accession,go_obo)]
b73v4_gramene_gaf[,db_object_id:=db_object_symbol]
b73v4_gramene_gaf[,db:="Gramene"]
b73v4_gramene_gaf[,qualifier:=""]
b73v4_gramene_gaf[,db_reference:=""]
b73v4_gramene_gaf[,with:=""]
b73v4_gramene_gaf[,db_object_name:=""]
b73v4_gramene_gaf[,db_object_synonym:=""]
b73v4_gramene_gaf[,db_object_type:="protein"]
b73v4_gramene_gaf[,taxon:="taxon:4577"]
b73v4_gramene_gaf[,date:="20200606"]
b73v4_gramene_gaf[,assigned_by:="Gramene"]
b73v4_gramene_gaf[,c(setdiff(gaf_cols,colnames(b73v4_gramene_gaf))):=""]

cat("removing redundancy\n")
b73v4_gramene_gaf_nonred <- rm_gaf_red(b73v4_gramene_gaf[aspect!="N/A"],go_obo)
dim(b73v4_gramene_gaf_nonred)
write_gaf(b73v4_gramene_gaf_nonred,"data/predicted/maize.B73.AGPv4.gramene.gaf")
