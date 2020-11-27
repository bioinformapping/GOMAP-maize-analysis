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

go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")
config

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


venn_term_acc = split(predData_dt[,.(inbred,aspect,source,term_accession)],by=c("inbred","aspect","source"),flatten=FALSE,keep.by = F,sorted = T)
plotCols=c("#fdd9d6","#d8d3dc","#d9e6f1")
fillScale=plotCols[c(1,3)]
b73v4 = lapply(venn_term_acc$B73v4,function(aspectVennData){
  data <- lapply(aspectVennData,function(x){
    unique(unlist(go_obo$ancestors[unlist(x)]))
  })
  venn.gTree <- getVenn(data)
  venn.gTree
})
ph207 = lapply(venn_term_acc$PH207,function(aspectVennData){
  data <- lapply(aspectVennData,function(x){
    unique(unlist(go_obo$ancestors[unlist(x)]))
  })
  venn.gTree <- getVenn(data)
  venn.gTree
})
#plotMultiVenn(plotFile = "paper/figures/cVsgTermVenn.png",b73Venns = b73v4,ph207Venns = ph207,legendCols = plotCols,mainTitle = "Annotated GO Terms")


distinctB73TermList = lapply(venn_term_acc$B73v4,function(x){
  allCommunityTerms = unique(unlist(get_ancestors(go_obo,x$Community$term_accession)))
  allGomapTerms = unique(unlist(get_ancestors(go_obo,x$GOMAP$term_accession)))
  uniqCommunityTerms = setdiff(allCommunityTerms,allGomapTerms)
  communityPredDt <- predData_dt[inbred=="B73v4" & source=="Community" & term_accession %in% uniqCommunityTerms]
  uniqGomapTerms = setdiff(allGomapTerms,allCommunityTerms)
  gomapPredDt <- predData_dt[inbred=="B73v4" & source=="GOMAP" & term_accession %in% uniqGomapTerms]
  rbind(communityPredDt,gomapPredDt)
})
distinctB73Terms= unique(rbindlist(distinctB73TermList)[,.(inbred,source,aspect,term_accession)])
distinctB73Terms[,specificity:=unlist(get_specificity(term_accession,go_obo))]
distinctB73Terms
distinctB73Terms
specBreaks = c(0,6,11,21,41,Inf)
specLbls = c("0-5","6-10","11-20","21-40",">40")
distinctB73Terms[,specClass:=factor(cut(specificity,breaks = specBreaks,labels = specLbls,include.lowest = T))]

#distinctSpecB73 <- distinctB73Terms[,list(count=ifelse(source=="GOMAP",as.integer(.N),as.integer(.N*-1))),by=list(inbred,source,aspect,specClass)]
distinctSpecB73 = distinctB73Terms[,list(count=.N),by=list(source,aspect,specClass)]
distinctSpecB73

uniqTermSpecPlot <- ggplot(distinctSpecB73,aes(x=specClass,y=count,color=source,shape=source)) + 
  geom_line(aes(group=source),size=1) + geom_point(size=3)  +
  facet_grid(.~aspect,labeller = labeller(aspect=aspect_lbl)) + 
  scale_color_brewer(type = "qual",palette = 4) +
  theme_linedraw() + theme_plot + theme(legend.position = "bottom") +
  labs(Y="Specificity",y="# of Terms",shape="Annotation\nSource",color="Annotation\nSource",title="Specificity of the GO terms that are unique to GOMAP or the Community in B73 Annotations")
uniqTermSpecPlot
ggsave("paper/figures/cVsgTermSpecificity.png",plot = uniqTermSpecPlot)
