source("R/plotTools.R")
source("R/gafTools.R")
source("R/oboTools.R")
library("data.table")
library("naturalsort")
library("kableExtra")
library("gtable")
library("yaml")

allMaizeSeqs <- fread("data/fasta/filt_fasta/faSeqData.tsv")
allMaizeSeqsSummary <-  allMaizeSeqs[,list(total=.N),by=inbred][order(inbred)]
allMaizeSeqsSummary
totalGenes <- allMaizeSeqsSummary$total
names(totalGenes) <- allMaizeSeqsSummary$inbred

go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yaml")

predData_list <- lapply(config$predicted,function(x){
  inputGaf <- file.path("data/clean",x$file)
  inbred = x$inbred
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:=x$source]
  input_gaf_data
})

predData_dt <- rbindlist(predData_list)
predData_dt[,source:=factor(source,levels = c("GOMAP","Community"),ordered = T)]
predData_dt[,list(count=.N),by=list(inbred,source)]
system.time({
  specList = get_specificity(unique(predData_dt$term_accession),go_obo)
  names(specList) <- unique(predData_dt$term_accession)
  print(length(unlist(specList[predData_dt$term_accession])))
  NROW(predData_dt)
  predData_dt[,specificity:=unlist(specList[term_accession])]
  
})


#dim(gafData_dt)
fwrite(predData_dt[,.(inbred,source,aspect,db_object_symbol,term_accession,specificity)],file = "data/metrics/gafSpecificity.tsv")

predDataSummary <- predData_dt[,list(annotated=length(unique(db_object_symbol)),numAnnots=.N,meanSpec=round(mean(specificity),2)),by=list(inbred,source)]
predDataSummary[,totGenes:=totalGenes[inbred]]
predDataSummary[,coverage:=round(annotated/totGenes*100,2)]
predDataSummary[,annotsPerGene:=round(numAnnots/annotated,2)]

setorderv(predDataSummary,c("inbred","source"))
setcolorder(predDataSummary,c("inbred","totGenes","annotated","coverage","numAnnots","annotsPerGene","meanSpec"))

curateData_list <- lapply(config$curated,function(x){
  inputGaf <- file.path("data/curated/",x$file)
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data
})

curatedData_dt <- rbindlist(curateData_list)
system.time({
  specList = get_specificity(unique(curatedData_dt$term_accession),go_obo)
  names(specList) <- unique(curatedData_dt$term_accession)
  specList[curatedData_dt$term_accession]
  curatedData_dt[,specificity:=unlist(specList[term_accession])]
})

curatedDataSummary <- curatedData_dt[,list(annotated=length(unique(db_object_symbol)),
                                          numAnnots=.N,
                                          meanSpec=round(mean(specificity),2)),
                                     by=list(inbred)]
curatedDataSummary[,totGenes:=totalGenes[inbred]]
curatedDataSummary[,coverage:=round(annotated/totalGenes*100,2)]
curatedDataSummary[,annotsPerGene:=round(numAnnots/annotated,2)]

gafDataSummary <- merge.data.table(predDataSummary,curatedDataSummary,by=c("inbred","totGenes"),suffixes = c(".pred",".curate"))
colnames(gafDataSummary)
summaryCols <- c("inbred","source","totGenes",naturalsort(colnames(gafDataSummary)[grep("coverage|annots|meanSpec",colnames(gafDataSummary))]))
cat(kable(gafDataSummary[,summaryCols,with=F],"markdown",row.names = F,format.args = list(big.mark=","),caption="Geneal Analysis Metrics"),file = "tables/README.md",sep = "\n",append=TRUE)

sourceColor = c("black","black")
fillScale = scale_fill_brewer(type="qual",palette = 4)
colorScale = scale_color_manual(values = sourceColor)
shapeSclae = scale_shape(solid = TRUE)
pos=position_dodge2(width = 0.9,preserve = "single")
pos2=position_dodge2(width = 0.5,preserve = "single",padding = 0.7)
pointSize=3
shapeYpos=-0.05

predCovSummary <- predData_dt[,list(annotated=length(unique(db_object_symbol))),by=list(inbred,source,aspect)]
predCovSummary[,totGenes:=totalGenes[inbred]]
predCovSummary[,coverage:=round(annotated/totGenes*100,2)]
predCovSummary[inbred=="B73v4"]

covPlot <- ggplot(predCovSummary, aes(x=inbred)) +
  geom_bar(mapping = aes(y=coverage,fill=inbred,color=source),stat = "identity",position = pos,size=0) +
  geom_point(mapping=aes(y=shapeYpos*max(coverage),shape=source),position = pos,size=pointSize) +
  fillScale + colorScale + shapeSclae +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) + 
  theme_linedraw() + theme_strip_nxaxis + 
  theme(legend.position = "none") +
  labs(y="Coverage (%)",x="Inbred Line")
#covPlot

annotsPerGene <- predData_dt[,list(count=.N),by=list(inbred,source,db_object_symbol,aspect)]
numAnnotData = annotsPerGene[,list(avg=mean(count),lower=t.test(count)$`conf.int`[1],upper=t.test(count)$`conf.int`[2]),by=list(inbred,source,aspect)]

numAnnotPlot <- ggplot(numAnnotData, aes(x=inbred,fill=inbred)) +
  geom_bar(aes(y=avg,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=shapeYpos*max(avg),shape=source),position = pos,size=pointSize) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) + 
  fillScale + colorScale + shapeSclae +
  theme_linedraw() + theme_nstrip_nxaxis + theme(legend.position = "bottom") +
  theme(legend.position = "none") +
  labs(y="# Annotations/Gene",x="Inbred Line",shape="Annotation\nSource",fill="Inbred\nLine")

specHistPlotData = predData_dt[,list(avg=mean(specificity),lower=t.test(specificity)$`conf.int`[1],upper=t.test(specificity)$`conf.int`[2]),by=list(inbred,source,aspect)]
#specHistPlotData = predData_dt[,list(avg=mean(specificity),lower=mean(specificity)-sd(specificity),upper=mean(specificity)+sd(specificity)),by=list(inbred,source,aspect)]
#specHistPlotData

specHistPlot <- ggplot(specHistPlotData, aes(x=inbred)) +
  geom_bar(aes(y=avg,fill=inbred,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=shapeYpos*max(avg),shape=source),position = pos,size=pointSize) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) + 
  fillScale + colorScale + shapeSclae +
  theme_linedraw() + theme_nstrip_xaxis +
  guides(color=F,fill=guide_legend(nrow = 2),shape=guide_legend(nrow = 2)) +
  theme(legend.position = "bottom") +
  labs(y="Specificity",x="Inbred Line",shape="Annotation\nSource",fill="Inbred\nLine")

#specHistPlot
#metricGrob = arrangeGrob(covPlot,numAnnotPlot,specHistPlot,heights = c(1.1,1,1.2),ncol = 1)

out =  lapply(list(covPlot,numAnnotPlot,specHistPlot),ggplotGrob)
g = do.call(rbind,out)
#plotGrob <- arrangeGrob(g,legends,heights = c(6,1),ncol = 1,padding = unit(2,"line"))
ggsave("figures/analMetrics.png",plot = g,width=7.5,height=6,units = "in")

