#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-22,2018-07-09_1,2018-07-09_2,2018-07-25,2018-10-16
# for : hot spot　recovery
# 参数说明
# Options
# -d/--dir        热点位点文件所在文件夹 　　　　
# -v/--vcf        vcf文件 　　　　　
# -o/--outdir     输出文件
options(warn = -1)
rm(list=ls()) 
package_list <- c("optparse","readxl","magrittr","stringr","tidyr","Biostrings","foreach")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("https://bioconductor.org/biocLite.R")
      biocLite(p)
    }
  }
}

sapply(package_list,FUN = function(n) require(n,character.only = T,quietly = T))

if(T){
  option_list <- list(make_option(c("-d", "--dir"), type="character",help = "mutation_hotspot Excel file dir"), 
                      make_option(c("-v", "--vcf"), type="character"), 
                      make_option(c("-o","--outdir"), type="character"))
  opts <- parse_args(OptionParser(option_list=option_list))
}


if(F){
      opts= list()
      opts$vcf = "test/Ct18100066_1121_10027805.snvIndel.raw.vcf"
      opts$dir = "/home/wqh/work/ctDNA_output_filter"
      opts$outdir = "."
}


s0 <- Sys.time() 
print(opts$vcf)
trans_REF_ALT <- function(data){
  Length = sapply(data[4:5], str_length)
  if(dim(data)[1]>1){
    Max_num =  apply(Length,1,which.max)
    Min_num =  apply(Length,1,which.min)
  }else{
    Max_num = which.max(Length)
    Min_num = which.min(Length)
  }
  data$pos = data$pos + 1
  for(i in 1:nrow(data)){
    tmp_data = data[i,]
    data[i,4:5][Max_num[i]] = sub(pattern = data[i,4:5][Min_num[i]],x = data[i,4:5][Max_num[i]],replacement = "")
    data[i,4:5][Min_num[i]] = "-"
  }
  return(data)
} 

filename = basename(opts$vcf) %>% str_split(pattern = "\\.",simplify = T) %>%'['(1)
col_names = "#chr	pos_start	pos_end	ref	alt	Gene	Type	CDSChange	exon	sequence	AA	X__1	chr	pos	ref.1	alt.1	d" %>% str_split(pattern = "\\s",simplify = T)

file = read_excel(sprintf("%s/mutation_hotspot.xlsx",opts$dir),sheet = 2)
file$`#chr` %>% str_replace(pattern="chr",replacement ="")  %>% as.numeric() -> file$`#chr`
file <- file[order(file$`#chr`,file$pos_start),]
file %<>% unite(col = "IID",c("#chr","pos_start"),sep = ":",remove = F)

vcf <- read.table(opts$vcf,stringsAsFactors = F)
colnames(vcf) <- c("chr","pos","a","ref","alt","b","c","d","e","f")

if(any(grepl(vcf$d,pattern = "INDEL"))) vcf[sapply(vcf$d, grepl, pattern="INDEL",simplify = T),] %<>% trans_REF_ALT()
vcf %<>% unite(col = "IID",c("chr","pos"),sep = ":",remove = F)
sub_vcf  = vcf[vcf$IID %in% file$IID, ]


if(nrow(sub_vcf)!=0){
  a <- function(i,sub_vcf,file){
    tmp_ref_alt = sub_vcf[i,5:6] %>% paste(collapse = ":")
    tmp_ref_alt_comlement = sub_vcf[i,5:6] %>% sapply(function(n) n %>% DNAString() %>% complement() %>% as.character) %>% paste(collapse = ":")
    tmp_file = subset(file,IID == sub_vcf[i,1]) %>% unite(col = "IID",c("ref","alt"),sep = ":",remove = F)
    tmp_hotspot = tmp_file[(tmp_file$IID %in% tmp_ref_alt)|(tmp_file$IID %in% tmp_ref_alt_comlement),-5]
    if (nrow(tmp_hotspot)==1) cbind(tmp_hotspot,sub_vcf[i,c(2,3,5,6,9)])}
  out = foreach(i=1:nrow(sub_vcf),.combine = rbind,sub_vcf,file) %do% a(i,sub_vcf = sub_vcf ,file = file)
  dat_out=NULL
  if(!is.null(out)) dat_out = out
  if(!is.null(dat_out)){
    write.table(col_names,paste0(opts$outdir,"/",filename,".hotspot.xls"),sep = "\t",quote = F,row.names = F,col.names = F)
    write.table(dat_out,paste0(opts$outdir,"/",filename,".hotspot.xls"),sep = "\t",col.names = F,append = T,quote = F,row.names = F)
  }else{
    write.table(col_names,paste0(opts$outdir,"/",filename,".hotspot.xls"),sep = "\t",quote = F,row.names = F,col.names = F)
  }
}else{
  write.table(col_names,paste0(opts$outdir,"/",filename,".hotspot.xls"),sep = "\t",quote = F,row.names = F,col.names = F)
}
 
Sys.time() - s0

