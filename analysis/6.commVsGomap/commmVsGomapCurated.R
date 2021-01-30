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
config <- read_yaml("config.yaml")
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

# Filtering predicted data for the gene ids that have annotations in curated data
curateData_dt = rbindlist(curateData_list)
curateFilter = unique(curateData_dt[,.(db_object_id,aspect)])
predDataFilt = merge(predData_dt,curateFilter)

# Combine the predicted and curated datasets
allAnnotData = rbind(predDataFilt,curateData_dt)

# Obtain all the terms till the root term for each GO term in the predicted and curated annotations
allAnnots = unique(allAnnotData[,.(term_accession)])
annotTermsAncest = allAnnots[,list(full_term_accession=go_obo$ancestors[[term_accession]]),by=term_accession]
system.time({
  all_term_accession <- unique(merge(allAnnotData,annotTermsAncest,allow.cartesian = T))
  all_term_accession
})


# Filter based on whether term is present in the curated data or not. We only 
# care about the terms that are in the curated annotations and no other terms 
curatedFilt = unique(all_term_accession[source=="Curated",.(db_object_symbol,aspect,full_term_accession,inbred)])
all_term_data = merge.data.table(all_term_accession,curatedFilt,by=c("db_object_symbol","aspect","full_term_accession","inbred"))
all_term_data[,annotPair:=paste(db_object_symbol,full_term_accession,sep="-")]
all_term_data[,inData:=T]

all_term_data


getAnnotCat = function(annotSource){
  gomap="GOMAP" %in% annotSource
  community= "Community" %in% annotSource
  curated="Curated" %in% annotSource
  if(gomap & community){
    return("Both")
  }else if(gomap & !community){
    return("GOMAP")
  }else if(!gomap & community){
    return("Community")
  }else{
    return("Curated")
  }
  return("NULL")
}

all_term_data

geneCols=c("db_object_symbol","aspect","inbred","source")
geneCounts = unique(all_term_data[,geneCols,with=F])[,list(annotClass=getAnnotCat(source)),by=list(db_object_symbol,aspect,inbred)][,list(variable="Genes",count=.N),by=list(aspect,inbred,annotClass)][order(aspect,inbred,annotClass)]
termCols=c("full_term_accession","aspect","inbred","source")
termCounts = unique(all_term_data[,termCols,with=F])[,list(annotClass=getAnnotCat(source)),by=list(full_term_accession,aspect,inbred)][,list(variable="GO Terms",count=.N),by=list(aspect,inbred,annotClass)][order(aspect,inbred,annotClass)]
selCols=c("db_object_symbol","aspect","full_term_accession","inbred","source")
annotCounts = unique(all_term_data[,selCols,with=F])[,list(annotClass=getAnnotCat(source)),by=list(db_object_symbol,full_term_accession,aspect,inbred)][,list(variable="Annotations",count=.N),by=list(aspect,inbred,annotClass)][order(aspect,inbred,annotClass)]

plotData = rbindlist(list(geneCounts,termCounts,annotCounts),use.names = T)
plotData[,prop:=count/sum(count),by=list(aspect,inbred,variable)]
plotData[,annotClass:=factor(annotClass,levels = c("Both","GOMAP","Community","Curated"))]
plotData[,variable:=factor(variable,levels = c("Genes","GO Terms","Annotations"))]
plotData

cols = brewer_pal(type = "qual",palette = 4)

p = ggplot(plotData,aes(x=inbred,y=prop,fill=annotClass)) +
  geom_bar(stat="identity") +
  facet_grid(variable~aspect,scales = "free",labeller = labeller(aspect=aspect_lbl)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cols(4)) +
  theme_linedraw(base_size = 12) + theme_plot + theme(legend.position = "bottom") +
  labs(x="Inbred Line",y="Proportion (%)",fill="Annotated By")
plot(p) 

ggsave("figures/gomapVsCommVsCurate.png",plot = p,width = 6,height = 7.5,units = "in",dpi = 300)

?muted
countTable = dcast.data.table(plotData,inbred+annotClass~aspect+variable,value.var = "count",fill = 0)#[order()]
cat(colSums(countTable[inbred=="B73v4",-1:-2,with=F]),sep=" & ")
cat(colSums(countTable[inbred!="B73v4",-1:-2,with=F]),sep=" & ")

kable(countTable,format = "latex",format.args = list(big.mark=","))
