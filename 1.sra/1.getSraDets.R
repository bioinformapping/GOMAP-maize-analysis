# Importing Libraries -----------------------------------------------------
library("data.table")
library("httr")
library("jsonlite")
library("xml2")
library("parallel")


queryUrl="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=((%22pacbio%20smrt%22[Platform])+OR+%22oxford%20nanopore%22[Platform])+AND+Embryophyta[Organism]+AND+wgs[strategy]&retmax=100000&format=json"
GET(queryUrl,write_disk("data/sra/sra_query.json",overwrite = F))
id_list = read_json("data/sra/sra_query.json")$esearchresult$idlist
names(id_list) = rep("id",length(id_list))

# Obtaining SRA Results ---------------------------------------------------
efetchUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
idlist1 = id_list
idlist1[["db"]] ="sra"
out = POST(url = efetchUrl,body = idlist1,encode = "form",write_disk("data/sra/sra_query.xml"), verbose())


# Process SRA XML ---------------------------------------------------------

tmpXmlDoc <- read_xml("data/sra/sra_query.xml")
expPkgs <- xml_find_all(tmpXmlDoc,".//EXPERIMENT_PACKAGE")
length(expPkgs)
out1 <- mclapply(expPkgs,function(y){
  experiment <- xml_attr(xml_find_all(y,".//EXPERIMENT"),"accession")
  study <- xml_attr(xml_find_all(y,".//STUDY"),"accession")
  study_title <- xml_text(xml_find_all(y,".//STUDY/DESCRIPTOR/STUDY_TITLE"),"accession")
  sci_name <- xml_text(xml_find_all(y,".//SCIENTIFIC_NAME"))
  runSet <- xml_find_first(y,".//RUN_SET")
  organization <- xml_text(xml_find_all(y,".//Organization/Name"))
  runDets <- xml_find_first(runSet,".//RUN")
  selCols <- c("accession"," published")
  accession <- xml_attr(runDets,"accession")
  pub_date <- xml_attr(runDets,"published")
  data.table(cbind(organization,experiment,study,study_title,sci_name,accession,pub_date))
},mc.cores = 8)
outDt <- rbindlist(out1,fill = T)
fwrite(outDt,file = "data/sra/sra_query.tsv",sep = "\t",row.names = F,quote = F)

