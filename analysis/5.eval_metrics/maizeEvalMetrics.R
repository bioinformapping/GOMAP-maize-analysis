library("data.table")
library("naturalsort")
library("xtable")
source("analysis/R/cafaEval.R")
source("analysis/R/aigoEval.R")
source("analysis/R/plotTools.R")
library("ggplot2")
library("gridExtra")


allMaizeSeqs <- fread("data/2.filt_fasta/faSeqData.tsv")
allMaizeSeqsSummary <-  allMaizeSeqs[,list(total=.N),by=inbred][order(inbred)]
allMaizeSeqsSummary
totalGenes <- allMaizeSeqsSummary$total
names(totalGenes) <- allMaizeSeqsSummary$inbred
totalGenes <- c(totalGenes,c("B73v3-Gold"=as.numeric(totalGenes["B73v3"])))
totalGenes

#B73v4 <- read_gaf("predicted/maize.B73.AGPv4.aggregate.gaf")

go_obo <- check_obo_data("data/go/go.obo")

inbredLines <- c("B73v3","B73v4","Mo17","PH207","W22")
predFiles <- dir("data/6.predicted/","aggregate",full.names = T)
names(predFiles) <- inbredLines
goldFiles <- dir("data/5.curated/",pattern = "curated",full.names = T)
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
#save(gafData_dt,file = "Rdata/gafCombined.data")
saveRDS(gafData_dt,"gafCombined.data")
gafData_dt

curatedGenes <- unique(gafData_dt[annotType=="curated"]$db_object_symbol)
gafEvalData <- gafData_dt[db_object_symbol %in% curatedGenes]
gafEvalData

system.time({
  gafCafaMetrics <- gafEvalData[,getCafaEval(.SD),by=list(inbred,db_object_symbol,aspect)]
})
gafCafaMetrics[,fscore:=2*(pr*rc/(pr+rc))]
gafCafaMetrics[is.nan(fscore),fscore:=0]
#gafCafaMetrics[num_gold>0 & num_pred==0]
gafCafaMetrics
    
curPredCount_dt <- gafCafaMetrics[num_gold>0 & num_pred>0,list(predCount=.N),by=list(inbred,aspect)][order(inbred)]
curCount_dt <- gafData_dt[annotType=="curated",list(totCount=length(unique(db_object_symbol))),by=list(inbred,aspect)]

curPredCount <- merge(curCount_dt,curPredCount_dt)
curPredCount[,curCount:=totCount-predCount]
curPredCount
curPredCount_melt <- melt(curPredCount,id.vars = c("inbred","aspect"))
curPredCount_melt

curPred_plot <- ggplot(curPredCount_melt[variable!="totCount"],aes(x=inbred,y=value,fill=variable)) + 
  geom_bar(stat="identity") + 
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  scale_fill_brewer(type="qual",palette = 2,labels=geneType_list) + 
  theme_linedraw(base_size = 12) + theme_plot + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  labs(y="# of Genes with Curated Terms",fill="Annotation Type")
                
curPred_plot

gafCafaMetrics[num_gold>0][,list(pr=mean(pr),rc=mean(rc)),by=list(inbred,aspect)]
gafCafaMetrics[,fscore:=2*(pr*rc/(pr+rc))]
gafCafaMetrics[is.nan(fscore),fscore:=0]
gafCafaMetrics

gafCafaMetrics[num_gold>0,list(fscore=2*(mean(pr)*mean(rc)/(mean(pr)+mean(rc)))),by=list(inbred,aspect)]

prRc_melt <- melt(gafCafaMetrics[num_gold>0 & num_pred>0],id.vars = c("inbred","db_object_symbol","aspect","evalType"),measure.vars = c("pr","rc"))
prRc_melt
prRc_melt
gafCafaMetrics


pr_plot <- ggplot(gafCafaMetrics[num_gold>0],aes(x=inbred,y=pr,fill=inbred)) + 
  geom_boxplot() + 
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  scale_fill_brewer(type="qual",palette = 1) +
  scale_color_brewer(type="qual",palette = 1) + 
  theme_linedraw(base_size = 12) + theme_plot + guides(fill=F) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) +
  labs(y="Precision",fill="Annotation Type")
                

rc_plot <- ggplot(gafCafaMetrics[num_gold>0],aes(x=inbred,y=rc,fill=inbred)) + 
  geom_boxplot() + 
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  scale_fill_brewer(type="qual",palette = 1) +
  scale_color_brewer(type="qual",palette = 1) + 
  theme_linedraw(base_size = 12) + theme_plot + guides(fill=F) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_blank(),strip.text.x = element_blank(),axis.title.x = element_blank()) +
  labs(y="Recall",fill="Annotation Type")
#prRc_plot

cafgaFscore <- gafCafaMetrics[num_gold>0,list(fscore=2*(mean(pr)*mean(rc)/(mean(pr)+mean(rc)))),by=list(inbred,aspect)] 
fscore_plot <- ggplot(cafgaFscore,aes(x=inbred,y=fscore,fill=inbred)) + 
  geom_bar(stat="identity") + 
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  scale_fill_brewer(type="qual",palette = 1) +
  scale_color_brewer(type="qual",palette = 1) + 
  theme_linedraw(base_size = 12) + theme_plot + guides(fill=F) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(y=bquote(F[max]))
                
#fscore_plot

plotGrob <- arrangeGrob(pr_plot,rc_plot,fscore_plot,ncol = 1,heights = c(1,1,1.2))
plot(plotGrob)

ggsave("paper/figures/cafaMetrics.png",plot = plotGrob,width=6.5,height=6,units = "in")

gafEvalData[db_object_symbol=="Zm00004b040638",getGamerEval(.SD),by=list(inbred,db_object_symbol,aspect)]

system.time({
  gafGamerMetrics <- gafEvalData[,getGamerEval(.SD),by=list(inbred,db_object_symbol,aspect)]
})

gafGamerMetrics[,fscore:=2*(pr*rc/(pr+rc))]
gafGamerMetrics[is.nan(fscore),fscore:=0]

gafCafaMetrics[num_gold>0,list(pr=mean(pr),rc=mean(rc),fscore=mean(fscore)),by=list(inbred,aspect)][inbred=="B73v3"]
gafGamerMetrics[num_gold>0,list(pr=mean(pr),rc=mean(rc),fscore=mean(fscore)),by=list(inbred,aspect)][inbred=="B73v3"]

2*(0.63*0.55)/(0.63+0.55)
gafGamerMetrics[db_object_symbol=="Zm00004b040638"]

gafEvalData[db_object_symbol=="Zm00004b040638",getAigoEval(.SD),by=list(inbred,db_object_symbol,aspect)]

system.time({
  gafAigoMetrics <- gafEvalData[,getAigoEval(.SD),by=list(inbred,db_object_symbol,aspect)]
})
gafAigoMetrics[,fscore:=2*(pr*rc/(pr+rc))]
gafAigoMetrics[is.nan(fscore),fscore:=0]

gafCafaMetrics[num_gold>0,list(pr=mean(pr),rc=mean(rc),fscore=mean(fscore)),by=list(inbred,aspect)][inbred=="B73v3"]
gafGamerMetrics[num_gold>0,list(pr=mean(pr),rc=mean(rc),fscore=mean(fscore)),by=list(inbred,aspect)][inbred=="B73v3"]
gafAigoMetrics[num_gold>0,list(pr=mean(pr),rc=mean(rc),fscore=mean(fscore)),by=list(inbred,aspect)][inbred=="B73v3"]
