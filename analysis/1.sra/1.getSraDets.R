library("data.table")
library("xml2")

xmlDoc <- xml2::read_xml("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=((%22pacbio%20smrt%22[Platform])+OR+%22oxford%20nanopore%22[Platform])+AND+Embryophyta[Organism]+AND+wgs[strategy]&retmax=100000")

datasetIds <- xml_text(xml_find_all(xmlDoc,".//Id"))
paste(datasetIds[1:10],collapse=",")
batchStarts <- seq(to = length(datasetIds),by = 100)
windWidth = 99
length(batchStarts)

runDetails <- lapply(batchStarts,function(x){
  startIdx=as.numeric(x)
  endIdx=startIdx+windWidth
  if(endIdx>length(datasetIds)){
    endIdx=length(datasetIds)
  }
  print(which(batchStarts == x))
  print(endIdx)
  outFile <- paste("analysis/1.sra/data/",startIdx,"-",endIdx,".tsv",sep="")
  if(!file.exists(outFile)){
    Sys.sleep(1)
    ids <- datasetIds[startIdx:endIdx]
    idTxt=paste(datasetIds[startIdx:endIdx],collapse = ",")
    tmpUrl <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=",idTxt,sep="")
    #print(tmpUrl)
    tmpXmlDoc <- read_xml(tmpUrl)
    expPkgs <- xml_find_all(tmpXmlDoc,".//EXPERIMENT_PACKAGE")
    out1 <- lapply(expPkgs,function(y){
      experiment <- xml_attr(xml_find_all(y,".//EXPERIMENT"),"accession")
      study <- xml_attr(xml_find_all(y,".//STUDY"),"accession")
      study_title <- xml_text(xml_find_all(y,".//STUDY/DESCRIPTOR/STUDY_TITLE"),"accession")
      sci_name <- xml_text(xml_find_all(y,".//SCIENTIFIC_NAME"))
      runSet <- xml_find_first(y,".//RUN_SET")
      organization <- xml_text(xml_find_all(y,".//Organization/Name"))
      #print(runSet)
      runDets <- xml_find_first(runSet,".//RUN")
      selCols <- c("accession"," published")
      accession <- xml_attr(runDets,"accession")
      pub_date <- xml_attr(runDets,"published")
      data.table(cbind(organization,experiment,study,study_title,sci_name,accession,pub_date))
    })
    outDt <- rbindlist(out1,fill = T)
    outDt[,id:=ids]
    fwrite(outDt,file = outFile,sep = "\t",row.names = F,quote = F)
    outFile
  }
})


