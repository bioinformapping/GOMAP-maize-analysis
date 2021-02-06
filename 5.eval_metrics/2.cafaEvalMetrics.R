source("R/cafaEval.R")
source("R/plotTools.R")
library("yaml")
library("data.table")
library("ggplot2")
library("naturalsort")
library("gtable")

allMaizeSeqs <- fread("data/fasta/filt_fasta/faSeqData.tsv")
allMaizeSeqsSummary <-  allMaizeSeqs[,list(total=.N),by=inbred][order(inbred)]
totalGenes <- allMaizeSeqsSummary$total
names(totalGenes) <- allMaizeSeqsSummary$inbred
totalGenes <- c(totalGenes,c("B73v3-Gold"=as.numeric(totalGenes["B73v3"])))

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

curateData_list <- lapply(config$curated,function(x){
  inputGaf <- file.path("data/curated/",x$file)
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:="curated"]
  input_gaf_data
})

curatedData_dt <- rbindlist(curateData_list)
curatedData_dt[,annotType:="curated"]
curUniq = unique(curatedData_dt[,.(inbred,db_object_symbol,aspect)])

prData_filt = merge(curUniq,predData_dt)

gafCafaPr <- prData_filt[,getCafaPr(term_accession,inbred,aspect,db_object_symbol,go_obo,curatedData_dt),by=list(source,inbred,aspect,db_object_symbol)]
gafCafaPr[,source:=factor(source,levels = c("GOMAP","Community"),ordered = T)]
fwrite(gafCafaPr,file = "data/metrics/gafCafaPr.tsv")

gafCafaRc_list <- apply(unique(prData_filt[,.(inbred,source)]),1,function(x){
  pred_data <- prData_filt[inbred==x["inbred"] & source==x["source"]]
  cur_data  <- curatedData_dt[inbred==x["inbred"]]
  gafCafaRc <- cur_data[,getCafaRc(term_accession,inbred,aspect,db_object_symbol,go_obo,pred_data),by=list(inbred,aspect,db_object_symbol)]
  gafCafaRc[,source:=x["source"]]
})

gafCafaRc <- rbindlist(gafCafaRc_list)
gafCafaRc[,source:=factor(source,levels = c("GOMAP","Community"),ordered = T)]
fwrite(gafCafaRc,file = "data/metrics/gafCafaRc.tsv")

sourceColor = c("black","black")
fillScale = scale_fill_brewer(type="qual",palette = 4)
colorScale = scale_color_manual(values = sourceColor)
shapeSclae = scale_shape(solid = TRUE)
pos=position_dodge2(width = 0.9,preserve = "single")
pos2=position_dodge2(width = 0.5,preserve = "single",padding = 0.7)
shapeYpos=-0.05

curGeneCount <- curatedData_dt[,list(curCount=length(unique(db_object_symbol))),by=list(inbred,aspect)]

pr_plotData <- gafCafaPr[,list(avg=mean(pr),
                                    lower=t.test(pr)$`conf.int`[1],
                                    upper=t.test(pr)$`conf.int`[2],
                                    count=length(unique(db_object_symbol))),
                              by=list(inbred,aspect,source)]

pr_plot <- ggplot(pr_plotData,aes(x=inbred,fill=inbred)) + 
  geom_bar(aes(y=avg,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=shapeYpos,shape=source),position = pos,size=3) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  fillScale + colorScale + shapeSclae +
  theme_linedraw() + theme_strip_nxaxis +
  theme(legend.position = "none") +
  labs(y="Precision",fill="Annotation Type")

rc_plotData <- gafCafaRc[,list(avg=mean(rc),
                                    lower=t.test(rc)$`conf.int`[1],
                                    upper=t.test(rc)$`conf.int`[2],
                                    count=sum(num_gold>0)),
                              by=list(inbred,aspect,source)]

rc_plot <- ggplot(rc_plotData,aes(x=inbred,fill=inbred)) + 
  geom_bar(aes(y=avg,color=source),stat = "identity",position = pos,size=0) +
  geom_errorbar(aes(ymin=lower,ymax=upper,color=source),position = pos2) +
  geom_point(aes(y=shapeYpos,shape=source),position = pos,size=3) +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl),scales="free") + 
  fillScale + colorScale + shapeSclae +
  theme_linedraw() + theme_nstrip_nxaxis + theme(legend.position = "none") +
  labs(y="Recall",fill="Annotation Type")

fMaxdata <- merge.data.table(pr_plotData[,.(inbred,source,aspect,avg)],
                             rc_plotData[,.(inbred,source,aspect,avg)],
                             by = c("inbred","source","aspect"),
                             suffixes = c(".pr",".rc"))

fMaxdata[,fmax:=2*(avg.pr*avg.rc)/(avg.pr+avg.rc)]

fMaxPlot <- ggplot(fMaxdata, aes(x=inbred)) +
  geom_bar(mapping = aes(y=fmax,fill=inbred,color=source),stat = "identity",position = pos,size=0) +
  geom_point(mapping=aes(y=shapeYpos,shape=source),position = pos,color="#000000",size=3) +
  fillScale + colorScale + shapeSclae +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) + 
  theme_linedraw() + theme_nstrip_xaxis + 
  theme(legend.position = "bottom") + 
  guides(color=F,fill=guide_legend(nrow = 2),shape=guide_legend(nrow = 2)) +
  labs(y=expression(italic(F[max])),x="Inbred Line",shape="Annotation\nSource",fill="Inbred\nLine")
fMaxPlot

out =  lapply(list(pr_plot,rc_plot,fMaxPlot),ggplotGrob)
g = do.call(rbind,out)
#plotGrob <- arrangeGrob(g,heights = c(6,1),ncol = 1,padding = unit(2,"line"))
ggsave("figures/cafaMetrics.png",plot = g,width=7.5,height=6,units = "in")