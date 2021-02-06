source("analysis/R/plotTools.R")
source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
source("analysis/R/goCompTools.R")
library("data.table")
library("naturalsort")
library("kableExtra")
library("gtable")
library("yaml")

## Loading the GO OBO data
go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")
config

## Loading the predicted data list
predData_list <- lapply(config$predicted[grep("B73v4",config$predicted)],function(x){
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

curateData_list <- lapply(config$curated[grep("B73v4",config$curated)],function(x){
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

totalCount = data.table(table(all_term_data$source))
totalCount

goldTerms = unique(curatedFilt[,.(aspect,full_term_accession)])


getUniq = function(full_term_accession_in,aspect_in,go_obo,totalCount){
  tmpDt = all_term_data[full_term_accession==unlist(full_term_accession_in) & aspect==unlist(aspect_in)]
  splitData = split(tmpDt[,.(source,db_object_symbol)],by=c("source"),flatten=F,keep.by = F)
  gomapCount = length(unique(splitData$GOMAP$db_object_symbol))
  gomapUniqCount = length(unique(setdiff(splitData$GOMAP$db_object_symbol,c(splitData$Community$db_object_symbol))))
  commCount = length(unique(splitData$Community$db_object_symbol))
  commUniqCount = length(unique(setdiff(splitData$Community$db_object_symbol,c(splitData$GOMAP$db_object_symbol))))
  curateCount = length(unique(splitData$Curated$db_object_symbol))
  curateUniqCount = length(unique(setdiff(splitData$Curated$db_object_symbol,c(splitData$GOMAP$db_object_symbol,splitData$Community$db_object_symbol))))
  out = unique(data.table(cbind(gomapCount,gomapUniqCount,commCount,commUniqCount,curateCount,curateUniqCount)))
  return(out)
}

uniqGeneCounts = goldTerms[,getUniq(full_term_accession,aspect,go_obo,totalCount),by=list(full_term_accession,aspect)]
uniqGeneCounts[full_term_accession=="GO:1901564"]
uniqGeneCounts

dim(uniqGeneCounts[curateUniqCount==0 & commUniqCount>0])
uniqGeneCounts[curateUniqCount==0 & gomapUniqCount>0,list(count=.N),by=list(aspect)]
uniqGeneCounts[curateUniqCount==0 & commUniqCount>0,list(count=.N),by=list(aspect)]
uniqGeneCounts[commUniqCount>0,list(count=.N),by=list(aspect)]
uniqGeneCounts[gomapUniqCount>0,list(count=.N),by=list(aspect)]

uniqGeneCounts[commCount+commUniqCount==0 & gomapUniqCount+gomapCount==0]
uniqGeneCounts[commCount+commUniqCount==0 & gomapUniqCount+gomapCount>0]
uniqGeneCounts[commCount+commUniqCount>0 & gomapUniqCount+gomapCount==0]


