library("data.table")
source("R/gafTools.R")
source("R/removeDuplicateAnnots.R")
source("R/getNonredDataset.R")
source("R/oboTools.R")

go_obo <- check_obo_data("data/go/go.obo")

## Clean the input dataset downloaded from maizeGDB

maizev4CrossGenes <-fread("data/maize/gene_model_xref_v4.txt")
selCols <- c("v4_gene_model","v3_gene_model","W22","PH207","Mo17b")
filtCrossGenes <- maizev4CrossGenes[,selCols,with=F]
colnames(filtCrossGenes) <- c("B73v4","B73v3","W22","PH207","Mo17")
fwrite(filtCrossGenes,"data/maize/maizeCrossRef.tsv",sep = "\t",row.names = F)

maizev4CrossGenes <- fread("data/maize/maizeCrossRef.tsv")

maizeCrossRef_melt <- melt.data.table(
  maizev4CrossGenes,
  id.vars = c("B73v3"),
  measure.vars = c("B73v4","Mo17","PH207","W22"),
  variable.name = "inbred",
  value.name = "geneModel"
)

splitComma <- function(data){
  geneModels <- unlist(strsplit(data,","))
  return(geneModels)
}

maizeCrossRef_genes <- maizeCrossRef_melt[,list(geneModel=splitComma(geneModel)),by=list(inbred,B73v3)]

curatedMaizeAnnot <- read_gaf("data/curated/maize_v3.gold.gaf")


write_gaf(curatedMaizeAnnot[grep("^GRM|^AC",db_object_symbol)],"data/curated/B73v3.curated.gaf")
annot_merge <- merge.data.table(maizeCrossRef_genes,curatedMaizeAnnot,by.x="B73v3",by.y="db_object_symbol")
annot_merge[,db_object_symbol:=geneModel]
annot_merge[,db_object_id:=geneModel]

out <- lapply(as.character(unique(annot_merge$inbred)),function(x){
    curGafFile <- paste("data/curated/",x,".curated.gaf",sep="")
    write_gaf(annot_merge[inbred==x,-c("B73v3","inbred","geneModel"),with=F],curGafFile)
})

