library("data.table")
library("naturalsort")
library("xtable")
source("R/cafaEval.R")
source("R/aigo_eval.r")
source("R/plot_tools.r")
library("ggplot2")
library("gridExtra")


allMaizeSeqs <- fread("maize/fa/faSeqData.tsv")
allMaizeSeqsSummary <-  allMaizeSeqs[,list(total=.N),by=inbred][order(inbred)]
allMaizeSeqsSummary
totalGenes <- allMaizeSeqsSummary$total
names(totalGenes) <- allMaizeSeqsSummary$inbred
totalGenes <- c(totalGenes,c("B73v3-Gold"=as.numeric(totalGenes["B73v3"])))
totalGenes

#B73v4 <- read_gaf("predicted/maize.B73.AGPv4.aggregate.gaf")

go_obo <- check_obo_data("go/go.obo")

inbredLines <- c("B73v3","B73v4","Mo17","PH207","W22")
predFiles <- dir("predicted","aggregate",full.names = T)
names(predFiles) <- inbredLines
goldFiles <- dir("curated",pattern = "curated",full.names = T)
names(goldFiles) <- inbredLines
goldFiles

gafData_list <- lapply(inbredLines,function(x){
  print(x)
  predGaf <- predFiles[[x]]
  print(predGaf)
  predGafData <- read_gaf(predGaf)
  predGafData[,inbred:=x]
  predGafData[,annotType:="predicted"]
  
  goldGaf <- goldFiles[[x]]
  print(goldGaf)
  goldGafData <- read_gaf(goldGaf)
  goldGafData[,inbred:=x]
  goldGafData[,annotType:="curated"]
  
  gafAll_dt <- rbindlist(list(predGafData,goldGafData),use.names = T)
})

gafData_dt <- rbindlist(gafData_list)
gafData_dt[,list(count=.N),by=list(inbred,annotType)]
save(gafData_dt,file = "Rdata/gafCombined.data")
gafData_dt

countPlotData <- gafData_dt[,list(geneCount=length(unique(db_object_symbol)),annotCount=.N),by=list(inbred,annotType)]
countPlotData
countPlotData[,annotsGene:=annotCount/geneCount]
#countPlotData[annotType=="predicted"]$geneCount = countPlotData[annotType=="predicted"]$geneCount - countPlotData[annotType=="curated"]$geneCount
countPlotData
countData_melt <- melt(countPlotData,id.vars = c("inbred","annotType"))
countData_melt

annotType_labs <- list(curated="Curated",predicted="Predicted")

geneCountPlot <- ggplot(countPlotData, aes(x=inbred,y=geneCount,fill=inbred)) +
  geom_bar(stat="identity",position = "dodge") +
  facet_grid(annotType~.,scales = "free",labeller = labeller(annotType=annotType_lbl)) +
  scale_color_brewer(type="qual",palette = 1) + 
  scale_fill_brewer(type="qual",palette = 1) + 
  scale_y_continuous(labels = comma) + 
  guides(fill=F) +
  theme_linedraw(base_size = 12) + theme_plot +
  labs(y="# of Annotated Genes",x="Inbred Line")
  
geneCountPlot

numAnnotPlot <- ggplot(countPlotData, aes(x=inbred,y=annotsGene,fill=inbred)) +
  geom_bar(stat="identity") +
  facet_grid(annotType~.,scales = "free",labeller = labeller(annotType=annotType_lbl)) +
  scale_fill_brewer(type="qual",palette = 1) +
  guides(fill=F) +
  theme_linedraw(base_size = 12) + theme_plot +
  labs(y="# of Annotations/Gene",x="Inbred Line")
numAnnotPlot

plotGrob <- arrangeGrob(geneCountPlot,numAnnotPlot,ncol = 2)
plot(plotGrob)

ggsave("/home/gokul/lab_data/papers/GOMAP-bioRxiv/figures/numGeneAnnot.png",plot = plotGrob,width=6.5,height=5,units = "in")
ggsave("plots/numGeneAnnot.png",plot = plotGrob,width=6.5,height=4,units = "in")
