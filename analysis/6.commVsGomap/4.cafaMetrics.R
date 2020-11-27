source("analysis/R/cafaEval.R")
source("analysis/R/plotTools.R")
library("yaml")
library("data.table")
library("ggplot2")
library("naturalsort")
library("gtable")
library("ontologySimilarity")
library("RColorBrewer")

go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")

commInbreds = unlist(lapply(config$predicted,function(x){
  if(x$source=="Community"){
    unique(x$inbred)
  }
}))
commInbreds

predData_list <- lapply(config$predicted,function(x){
  inputGaf <- file.path("data/clean",x$file)
  inbred = x$inbred
  if(inbred %in% commInbreds){
    input_gaf_data <- read_gaf(inputGaf)
    input_gaf_data[,inbred:=x$inbred]
    input_gaf_data[,source:=x$source]
    input_gaf_data
  }else{
    NULL
  }
})

predData_dt <- rbindlist(predData_list)
predData_dt
dcast(predData_dt[,list(count=length(unique(db_object_symbol))),list(inbred,aspect,source)],aspect~inbred+source)

curateData_list <- lapply(config$curated,function(x){
  inputGaf <- file.path("data/curated/",x$file)
  if(x$inbred %in% commInbreds){
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:="curated"]
  input_gaf_data
  }else{
    NULL
  }
})

curatedData_dt <- rbindlist(curateData_list)


curatedData_dt[,annotType:="curated"]
predData_dt[,annotType:="predicted"]

gafData_dt <- rbind(predData_dt,curatedData_dt)
gafData_dt[,list(count=.N),by=list(inbred,source,annotType)]

communityGenes <- unique(predData_dt[source=="Community",.(inbred,db_object_symbol,aspect)])
curatedGenes <- merge.data.table(communityGenes,unique(curatedData_dt[,.(inbred,db_object_symbol,aspect)]))
curatedGenes


gafEvalData[source=="GOMAP",list(count=.N),by=list(inbred,db_object_symbol,aspect)]

gafCafaMetrics_list <- mclapply(unique(predData_dt$source),function(x){
  print(x)
  gafPredData = merge.data.table(curatedGenes,gafData_dt[source==x])
  gafCurData  = merge.data.table(curatedData_dt,unique(gafPredData[,.(inbred,aspect)]))
  gafEvalData = rbind(gafCurData,gafPredData)
  print(gafEvalData[,list(count=.N),by=list(inbred,annotType)])
  #gafCafaMetrics <- gafEvalData[,getCafaEval(.SD),by=list(inbred,db_object_symbol,aspect)]
  gafCafaPr <- gafPredData[,getCafaPr(term_accession,inbred,aspect,db_object_symbol,go_obo,curatedData_dt),by=list(source,inbred,aspect,db_object_symbol)]
  gafCafaRc <- curatedData_dt[,getCafaRc(term_accession,inbred,aspect,db_object_symbol,go_obo,gafPredData),by=list(inbred,aspect,db_object_symbol)]
  print(gafCafaRc)
  print(gafCafaPr)
  gafCafaMetrics = merge.data.table(gafCafaRc,gafCafaPr,all = T)
  gafCafaMetrics[,source:=x]
  gafCafaMetrics
},mc.cores=4)
gafCafaMetrics = rbindlist(gafCafaMetrics_list)
gafCafaMetrics[num_gold>0 & num_pred>0,list(pr=mean(pr),rc=mean(rc)),by=list(inbred,aspect,source)]

dcast(gafCafaMetrics[,list(.N),by=list(inbred,aspect,source)],inbred+aspect~source)


pr_rc_melt = melt.data.table(gafCafaMetrics[num_gold>0 & num_pred>0],id.vars = c("inbred","source","aspect","db_object_symbol"),measure.vars = c("pr","rc"))

pr_rc_plotData = pr_rc_melt[,list(avg=mean(value),lower=t.test(value)$`conf.int`[1],upper=t.test(value)$`conf.int`[2]),by=list(inbred,source,aspect,variable)]
sourceColor = c("black","black")
fillScale = scale_fill_manual(values = brewer.pal(n = 5,name = "Pastel1")[c(2,4)])
colorScale = scale_color_manual(values = sourceColor)
pos=position_dodge2(width = 0.9,preserve = "single")
pos2=position_dodge2(width = 0.5,preserve = "single",padding = 0.7)

pr_plot <- ggplot(pr_rc_plotData[variable=="pr"],aes(x=inbred,fill=inbred)) +
  geom_bar(aes(y=avg,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=0,shape=source),position = pos,size=2) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") +
  fillScale + colorScale +
  theme_linedraw() + theme_strip_nxaxis +
  theme(legend.position = "none") +
  labs(y="Precision",fill="Annotation Type")

rc_plot <- ggplot(pr_rc_plotData[variable=="rc"],aes(x=inbred,fill=inbred)) +
  geom_bar(aes(y=avg,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=0,shape=source),position = pos,size=2) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") +
  fillScale + colorScale +
  theme_linedraw() + theme_nstrip_nxaxis + theme(legend.position = "none") +
  labs(y="Recall",fill="Annotation Type")

fMaxdata <- gafCafaMetrics[num_gold>0 & num_pred>0,list(fmax=2*(mean(pr)*mean(rc)/(mean(pr)+mean(rc)))),by=list(inbred,source,aspect)]
fMaxPlot <- ggplot(fMaxdata, aes(x=inbred)) +
  geom_bar(mapping = aes(y=fmax,fill=inbred,color=source),stat = "identity",position = pos,size=0) +
  geom_point(mapping=aes(y=0,shape=source),position = pos,color="#000000",size=2) +
  fillScale +
  colorScale +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) +
  theme_linedraw() + theme_nstrip_xaxis +
  theme(legend.position = "bottom") +
  guides(color=F,fill=guide_legend(nrow = 2)) +
  labs(y=expression(italic(F[max])),x="Inbred Line",shape="Annotation\nSource",fill="Inbred\nLine")
fMaxPlot

out =  lapply(list(pr_plot,rc_plot,fMaxPlot),ggplotGrob)
g = do.call(rbind,out)
ggsave("paper/figures/communityMetrics.png",plot = g,width=7.5,height=6,units = "in")
