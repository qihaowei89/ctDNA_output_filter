#!/usr/bin/env Rscript
# by: qihao
# date : 2018-09-30,2019-10-10
# for :  filter and reorder annotate.merge.filter.xls file


options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr","tidyr","Biostrings")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    if(!suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("http://bioconductor.org/biocLite.R")
      biocLite(p)
    }
  }
}

#--------------------------------------------- load required packages ---------------------------------------------#
sapply(package_list,FUN = function(n) require(n,character.only = T,quietly = T))

#--------------------------------------------- sub function ---------------------------------------------#
trans_3to1 <- function(n){ #蛋白质三字母到单字母
  n %<>% gsub(pattern = "p.",replacement = "")
  for (i in names(AMINO_ACID_CODE)) {
    n %<>%gsub(pattern =AMINO_ACID_CODE[i],replacement = i ) 
  }
  n %<>% gsub(pattern = "Ter",replacement = "*")
  n[grep(pattern = "%3D",n)] %<>% gsub(pattern = "([A-Z])([0-9]+)(.*)",replacement = "\\1\\2\\1")
  return(n)
}

#--------------------------------------------- get parameters ---------------------------------------------#
if(T){
  option_list <- list(make_option(c("-f", "--file"  ), type="character"),
                      make_option(c("-s", "--hotspot" ), type="character"),
                      make_option(c("-d","--dir"),type = "character"),
                      make_option(c("-o", "--outdir"), type="character"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "Usage: %prog [options]",add_help_option = T))
  Need_remove_site_xlsx = file.path(opts$dir,"Somtic_File_RemoveSite.xlsx")   
  Oncokb_annotate_xlsx  = file.path(opts$dir,"Oncokb-annotate.xlsx")  
  CKB_anV1_xlsx         = file.path(opts$dir,"CKB-anV1.xlsx") 
}

#--------------------------------------------- parameters for test ---------------------------------------------#
if(F){
  opts = list()
  opts$file    = "/home/wqh/workdir/ctDNA_output_filter/Ct1809_1113_NCCL1817_B.annotate.merge.filter.xls"
  opts$hotspot = "Ct1807_1113_HD_780_4.hotspot.xls"
  opts$dir = "RemoveSite_and_db"
  opts$outdir  = "test"
  Need_remove_site_xlsx = file.path(opts$dir,"Somtic_File_RemoveSite.xlsx")   
  Oncokb_annotate_xlsx  = file.path(opts$dir,"Oncokb-annotate.xlsx")  
  CKB_anV1_xlsx         = file.path(opts$dir,"CKB-anV1.xlsx") 
}


if(length(opts)!=5){
  cat("Usage:\n Rscript Somtic_RemoveSite.R/ \n    -f Ct1809_1113_NCCL1817_B.annotate.merge.filter.xls/ \n    -s Ct1809_1113_NCCL1817_B.hotspot.xls/ \n    -d RemoveSite_and_db/\n    -o test\n")
  cat(sprintf("3 options, only get %d\n",length(opts)-1))
  q()
}


#--------------------------------------------- 1.Oncokb and CKB db ---------------------------------------------#
s0 <- proc.time()
cat("step1: load Oncokb and CKB database\n")
Oncokb_annotate_db = read_xlsx(Oncokb_annotate_xlsx,range = cell_cols("D:I"))[,-2]  %>%  unite(col = "IID",c("Gene","Protein_Change"),sep = ":")
CKB_anV1_db = read_xlsx(CKB_anV1_xlsx,range = cell_cols("A:E"))[,-3] %>%  unite(col = "IID",c("Gene","Variant"),sep = ":")

# read raw file and order columns 
col_names       = file(opts$file,"r") %>% readLines(n = 1) %>% str_split(pattern = "\t",simplify = T)
file            = read.table(opts$file,header = F,sep = "\t",stringsAsFactors = F,col.names = col_names,check.names = F)
file$HotSpot    = rep("No",nrow(file))
file$Confidence = rep("No",nrow(file))
file$cosmic_id  = file$Cosmic82_coding %>% gsub(pattern = "^ID=",replacement = "")
file$cosmic_FATHMM_prediction       = rep(".",nrow(file))
file$cosmic_FATHMM_prediction_score = rep(".",nrow(file))
file_order_col  = file[,c(1,5,6,7,8,9,2,4,10,13,14,15,16,17:29,31:33,38,45:49)]


#--------------------------------------------- 2.Filter gene and regions ---------------------------------------------#
cat("step2: filter regions and target genes\n")
filter_regions = c("intronic", "UTR5", "UTR3", "ncRNA_intronic", "ncRNA_exonic", "intergenic")
str = "ABL1 /AKT1 /ALK /APC /AR /ARAF /ATM /BRAF /CCND1 /CDH1 /CDK4 /CDK6 /CDKN1A /CDKN2A /CTNNB1 /DDR2 /EGFR /ERBB2 /ERBB3 /ERBB4 /ESR1 /FBXW7 /FGFR1 /FGFR2 /FGFR3 /FGFR4 /FLT3 /GNA11 /GNAQ /GNAS /HRAS /IDH1 /IDH2 /JAK1 /JAK2 /JAK3 /KDR /KIT /KRAS /MAP2K1 /MAP2K2 /MAPK1 /MTOR /NF1 /NF2 /NRAS /NTRK1 /NTRK2 /NTRK3 /PDGFRA /PDGFRB /PIK3CA /POLE /PTCH1 /PTEN /RB1 /RET /ROS1 /SMAD4 /SMARCA4 /SMO /STAT3 /STK11  /TP53 /TSC1 /TSC2 /VHL"
filter_genes = sapply(str_split(str,pattern = "/"), function(n) trimws(n), simplify = T)

file_region_and_gene_filterd = file_order_col[((file_order_col$Gene %in% filter_genes) & !(file_order_col$Region %in% filter_regions)) | (file_order_col$Gene == "MET" & !(file_order_col$Region %in% filter_regions[-1])) | (file_order_col$Gene == "TERT" & (file_order_col$Start %in% c(1295228,1295250))), ]
file_region_and_gene_filterd$VEP_HGVSp %<>% trans_3to1()
file_region_and_gene_filterd %<>% unite(col = "IID",c("Gene","VEP_HGVSp"),sep = ":")
file_region_and_gene_filterd %<>%  merge(Oncokb_annotate_db,all.x=T) %>% merge(CKB_anV1_db,all.x=T) %>% separate(col = "IID",c("Gene","Protein_change"),sep = ":")
file_region_and_gene_filterd$Protein_change %<>% gsub(pattern = "([A-Z]+.*)",replacement = "p.\\1")

# reorder columns 
file_region_and_gene_filterd %<>% '['(c(3:8,32:34,1,9:11,2,12,13,16,14,15,31,30,36,35,37,38,39,17:29))
colnames(file_region_and_gene_filterd)[22:24] %<>% gsub(pattern = "(.*)",replacement = "Oncokb_\\1")
colnames(file_region_and_gene_filterd)[25:26] %<>% gsub(pattern = "(.*)",replacement = "CKB_\\1"   )

file_region_and_gene_filterd$Oncokb_Oncogenicity[is.na(file_region_and_gene_filterd$Oncokb_Oncogenicity)] = "."
file_region_and_gene_filterd$Oncokb_Mutation_Effect[is.na(file_region_and_gene_filterd$Oncokb_Mutation_Effect)] = "."
file_region_and_gene_filterd$Oncokb_PMIDs_for_Mutation_Effect[is.na(file_region_and_gene_filterd$Oncokb_PMIDs_for_Mutation_Effect)] = "."
file_region_and_gene_filterd$CKB_Protein_effect[is.na(file_region_and_gene_filterd$CKB_Protein_effect)] = "."
file_region_and_gene_filterd$CKB_Variant_description[is.na(file_region_and_gene_filterd$CKB_Variant_description)] = "."

# order output file 
file_out = file_region_and_gene_filterd %$% .[order(Chromosome,Start),] 


#---------------------------------------------  3.delete common sites ---------------------------------------------#
cat("step3: delete common sites\n")
delete_site  = read_xlsx(Need_remove_site_xlsx,sheet = 1,col_names = F) %>% unite(col="IID",c("X__1","X__2","X__3","X__4","X__5","X__6"),sep = ":")
file_out_tmp = file_out[,c(2:6,10)] %>% unite(col="IID",c(Chromosome,Start,End,Ref,Alt,Gene),sep = ":")
file_out = file_out[!(file_out_tmp$IID %in% delete_site$IID),]


# #删除常见突变,2018-07-18
# chrs  = c(4,4,4,5,5,5,5,5,6,6,7,8,9,9,9,9,9,9,10,11,11,11,12,14,19)
# Start = c(55972974,55561719,55561864,176520280,176520283,7870973,149510141,1295349,36651971,117647482,128829040,38285914,80537095,80537112,139399409,139399420,98209594,98211549,89685271,69465988,69466001,534295,56481662,105241411,1223125)
# End   = c(55972974,55561722,55561864,176520280,176520283,7870973,149510141,1295349,36651971,117647482,128829042,38285916,80537095,80537112,139399411,139399420,98209594,98211549,89685271,69465990,69466001,534297,56481662,105241411,1223125)
# Ref   = c("T","CCAT","A","C","G","A","T","A","C","T","GCT","TCA","G","T","CAC","C","G","G","T","GAG","A","CCA","G","A","C")
# Alt   = c("A","-",   "-","G","C","G","G","G","A","-",  "-",  "-","T","A",  "-","G","A","-","-",  "-","T",  "-","A","C","G")
# filter = data.frame(Chr=chrs,Start=Start,End=End,Ref=Ref,Alt=Alt,stringsAsFactors = F)
# for (i in 1:dim(filter)[1]) {
#   filter_row = filter[i,] 
#   delete = which((out$V1 %in% filter_row$Chr)& (out$V2 %in% filter_row$Start)&(out$V3 %in% filter_row$End)&(out$V4 %in% filter_row$Ref)&(out$V5 %in% filter_row$Alt))
#   if(length(delete)!=0) out = out[-delete,]
# }


#---------------------------------------------  4.search hotspot site   ---------------------------------------------#
cat("step4: search hotspot site\n")
hotspot_col_names   = file(opts$hotspot,"r") %>% readLines(n = 1) %>% str_split(pattern = "\t",simplify = T)
hotspot_file        = read.table(opts$hotspot,header = F,sep = "\t",stringsAsFactors = F,col.names = hotspot_col_names,check.names = F)[,-c(7:9,12,17)] %>% unite(col="IID",c("#chr","pos_start","pos_end","ref","alt","Gene","AA"),sep = ":")
file_out_tmp        = file_out[,c(2:6,10,13,14)] %>% unite(col = "IID",c("Chromosome","Start","End","Ref","Gene","Protein_change"),sep = ":")
check_if_in_hotspot = file_out_tmp$IID %in% hotspot_file$IID
file_out$HotSpot    = ifelse(check_if_in_hotspot,"Yes","No") 


#---------------------------------------------  5.write data to a file  ---------------------------------------------#
file_out_names = gsub(x=basename(opts$file),pattern="\\.xls",replacement=".final.xls")
cat(sprintf("step5: write data to : %s\n",file_out_names))
colnames(file_region_and_gene_filterd)[c(2,10:12)] = c("Chr","Gene_refGene","GeneDetail_ref","GeneExonicFunc_refGene")
write.table(x = file_out,file = file.path(opts$outdir,file_out_names),quote = F,sep = "\t",col.names = T,row.names = F)

s1 <- proc.time()
cat(sprintf("Complete !  Total Time: %.2fs\n",(s1-s0)[3]))

