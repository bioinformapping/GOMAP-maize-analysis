library("Biostrings")
library("kableExtra")
library("data.table")


inbredLines <- c("B73v3","B73v4","W22","PH207","Mo17")
faSeq_list <- lapply(inbredLines,function(x){
  print(x)
  inputFa <- paste("data/fasta/filt_fasta/",x,".fa",sep="")
  print(inputFa)
  input_seq <- Biostrings::readAAStringSet(inputFa)
  outDt <- data.table(cbind(faId=gsub(" .*","", names(input_seq)),seqLen=width(input_seq)    ))
  outDt[,inbred:=x]
  outDt
})
faSeq_dt <- rbindlist(faSeq_list)
faSeq_dt[,seqLen:=as.numeric(seqLen)]

faSeqSummary <- faSeq_dt[,list(numGenes=.N,totAALen=sum(seqLen),minAALen=min(seqLen),meanAAL=round(mean(seqLen),2),medAAL=median(seqLen),maxAALen=max(seqLen),smallGeneProp=round(sum(seqLen<50)/.N*100,2)),by=inbred]
setorderv(faSeqSummary,"inbred")
cat(kable(faSeqSummary,"markdown",row.names = F,format.args = list(big.mark=","),caption = "Fasta Sequence Summary"),file = "tables/README.md",sep = "\n",append=TRUE)
fwrite(faSeq_dt,file = "data/fasta/filt_fasta/faSeqData.tsv",sep = "\t",row.names = F)
