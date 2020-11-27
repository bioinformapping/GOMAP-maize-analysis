source("analysis/R/plotTools.R")
source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
source("analysis/R/goCompTools.R")
library("data.table")
library("naturalsort")
library("kableExtra")
library("gtable")
library("yaml")

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

## Loading the GO OBO data
go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")
config

## Loading the predicted data list
predData_list <- lapply(config$predicted[grep("B73v4|PH207",config$predicted)],function(x){
  inputGaf <- file.path("data/clean",x$file)
  inbred = x$inbred
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:=x$source]
  input_gaf_data
})

predData_dt <- rbindlist(predData_list )
unique(predData_dt[,.(inbred,source)])
predData_dt[,source:=factor(source,levels = c("GOMAP","Community"),ordered = T)]

curateData_list <- lapply(config$curated[grep("B73v4|PH207",config$curated)],function(x){
  inputGaf <- file.path("data/curated/",x$file)
  inbred = x$inbred
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:="Curated"]
  input_gaf_data
})

curateData_dt = rbindlist(curateData_list)
curateFilter = unique(curateData_dt[,.(db_object_id,aspect)])
predDataFilt = merge(predData_dt,curateFilter)

allAnnotData = rbind(predDataFilt,curateData_dt)

allAnnots = unique(allAnnotData[,.(term_accession)])
annotTermsAncest = allAnnots[,list(full_term_accession=go_obo$ancestors[[term_accession]]),by=term_accession]

system.time({
  all_term_accession <- unique(merge(allAnnotData,annotTermsAncest,allow.cartesian = T))
  all_term_accession
})

curatedFilt = unique(all_term_accession[source=="Curated",.(db_object_symbol,aspect,full_term_accession)])
curatedFilt

all_term_data = merge.data.table(all_term_accession,curatedFilt,by=c("db_object_symbol","aspect","full_term_accession"))
all_term_data[,annotPair:=paste(db_object_symbol,full_term_accession,sep="-")]


## Defining Fill Scales for pie chart
fillScale = brewer_pal(type = "qual",palette = 4)

## Splitting the data into hierarchical categories
vennData = split(all_term_data[,.(inbred,aspect,source,db_object_symbol)],by=c("inbred","aspect","source"),flatten=FALSE,keep.by = F,sorted = T)
## Generating the Venn data for B73 and PH207
#fillScale = alpha(c("red","yellow"),0.3)
plotCols=c("#66c2a5","#fc8d62","#8da0cb")
fillScale=plotCols
b73v4Genes = lapply(vennData$B73v4,function(aspectVennData){
  data <- lapply(aspectVennData,unlist)
  venn.gTree <- getVenn(data)
  venn.gTree
})

ph207Genes = lapply(vennData$PH207,function(aspectVennData){
  data <- lapply(aspectVennData,unlist)
  venn.gTree <- getVenn(data)
  venn.gTree
})

## Plotting the numbers for
plotMultiVenn(plotFile = "paper/figures/curateGeneVenn.png",
              b73Venns = b73v4Genes,
              ph207Venns = ph207Genes,
              legendCols = plotCols,
              mainTitle = "Curated Genes")

#system.time({
#  all_term_data <- predData_dt[sample(NROW(predData_dt),10000),
#    list(term_accession=unique(go_obo$ancestors[term_accession])),
#    by=list(inbred,aspect,db_object_symbol,source)]
#})


venn_term_acc = split(all_term_data[,.(inbred,aspect,source,full_term_accession)],by=c("inbred","aspect","source"),flatten=FALSE,keep.by = F,sorted = T)
#plotCols=c("#fdd9d6","#d8d3dc","#d9e6f1")
#fillScale=plotCols[c(1,3)]

b73v4Terms = lapply(venn_term_acc$B73v4,function(aspectVennData){
  data = lapply(aspectVennData,function(x){
    out1 = unique(unlist(x))
    out2 = unique(unlist(get_ancestors(go_obo,out1)))
    out2
  })
  venn.gTree <- getVenn(data)
  venn.gTree
})
ph207Terms = lapply(venn_term_acc$PH207,function(aspectVennData){
  data = lapply(aspectVennData,function(x){
    out1 = unique(unlist(x))
    out2 = unique(unlist(get_ancestors(go_obo,out1)))
    out2
  })
  venn.gTree <- getVenn(data)
  venn.gTree
})
plotMultiVenn(plotFile = "paper/figures/curateTermVenn.png",
              b73Venns = b73v4Terms,
              ph207Venns = ph207Terms,
              legendCols = plotCols,
              mainTitle = "Curated GO Terms")

venn_annot_pair = split(all_term_data[,.(inbred,aspect,source,annotPair)],by=c("inbred","aspect","source"),flatten=FALSE,keep.by = F,sorted = T)

b73v4Annots = lapply(venn_annot_pair$B73v4,function(aspectAnnotData){
  data = lapply(aspectAnnotData,unlist)
  venn.gTree <- getVenn(data)
  venn.gTree
})
ph207Annots = lapply(venn_annot_pair$PH207,function(aspectAnnotData){
  data = lapply(aspectAnnotData,unlist)
  venn.gTree <- getVenn(data)
  venn.gTree
})

plotInbredVenn(plotFile = "paper/figures/curatedB73v4Venn.png",geneVenns = b73v4Genes,termVenns = b73v4Terms,annotVenns = b73v4Annots,legendCols = plotCols,mainTitle = "B73v4 - GOMAP vs Community vs Curated")
plotInbredVenn(plotFile = "paper/figures/curatedPH207Venn.png",geneVenns = ph207Genes,termVenns = ph207Terms,annotVenns = ph207Annots,legendCols = plotCols,mainTitle = "PH207 - GOMAP vs Community vs Curated")

plotInbredVenn(plotFile = "../GOMAP-Plant-Methods/figures/curatedB73v4Venn.png",geneVenns = b73v4Genes,termVenns = b73v4Terms,annotVenns = b73v4Annots,legendCols = plotCols,mainTitle = "B73v4 - GOMAP vs Community vs Curated")
plotInbredVenn(plotFile = "../GOMAP-Plant-Methods/figures/curatedPH207Venn.png",geneVenns = ph207Genes,termVenns = ph207Terms,annotVenns = ph207Annots,legendCols = plotCols,mainTitle = "PH207 - GOMAP vs Community vs Curated")

