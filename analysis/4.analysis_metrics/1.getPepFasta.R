library("data.table")

if(!file.exists("data/fasta/raw_fasta/B73v3.fa.gz") & !file.exists("data/fasta/raw_fasta/B73v3.fa")){
  download.file("https://ftp.maizegdb.orgGDB/FTP/B73_RefGen_v3/Zea_mays.AGPv3.22.pep.all.fa.gz",destfile = "data/fasta/raw_fasta/B73v3.fa.gz")
}

### Downloading B73v4 Peptide fasta
if(!file.exists("data/fasta/raw_fasta/B73v4.fa.gz") & !file.exists("data/fasta/raw_fasta/B73v4.fa")){
  download.file("ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/zea_mays/pep/Zea_mays.AGPv4.pep.all.fa.gz",destfile = "data/fasta/raw_fasta/B73v4.fa.gz")
}

### Downloading PH207 Peptide fasta
if(!file.exists("data/fasta/raw_fasta/PH207.fa.gz") & !file.exists("data/fasta/raw_fasta/PH207.fa")){
  download.file("https://ftp.maizegdb.orgGDB/FTP/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0/Zm00008a.protein.fa.gz",destfile = "data/fasta/raw_fasta/PH207.fa.gz")
}

### Downloading Mo17 Peptide fasta
if(!file.exists("data/fasta/raw_fasta/Mo17.fa.gz") & !file.exists("data/fasta/raw_fasta/Mo17.fa")){
  download.file("https://ftp.maizegdb.orgGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/Zm00014a.proteins.fa.gz",destfile = "data/fasta/raw_fasta/Mo17.fa.gz")
}

### Downloading W22 Peptide fasta
if(!file.exists("data/fasta/raw_fasta/W22.fa.gz") & !file.exists("data/fasta/raw_fasta/W22.fa")){
  download.file("https://ftp.maizegdb.orgGDB/FTP/Zm-W22-REFERENCE-NRGENE-2.0/Zm00004b.protein.fa.gz",destfile = "data/fasta/raw_fasta/W22.fa.gz")
}