# 1 `src/hotspot.R`  
## 1.1 依赖:    

`optparse`, `readxl`, `magrittr`, `stringr`, `tidyr`, `Biostrings`, `foreach`包


## 1.2 用法:    
首次可运行 `Rscript hotsport.R -h `自动安装依赖的包      

显示如下结果表示依赖安装完成:   
```
    optparse     readxl   magrittr    stringr      tidyr Biostrings    foreach 
      TRUE       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
Usage: src/hotspot.R [options]


Options:
	-d DIR, --dir=DIR
		mutation_hotspot Excel file dir

	-v VCF, --vcf=VCF
		

	-o OUTDIR, --outdir=OUTDIR
		

	-h, --help
		Show this help message and exit
```

使用示例:   
`Rscript src/hotspot.R \
    -d ctDNA_output_filter/ \
    -v test/Ct18100047_1121_BB18100007.snvIndel.raw.vcf \
    -o ./ `

## 1.3 输入/输出文件   
- 输入: **`xxx.snvIndel.raw.vcf`**   
- 输出: **`xxx.hotspot.xls`**   

## 1.4 输出文件格式   
```
#chr	pos_start	pos_end	ref	alt	Gene	Type	CDSChange	exon	sequence	AA	X__1	chr	pos	ref.1	alt.1	d
7	55249071	55249071	C	T	EGFR	nonSynonymous_Substitution	NM_005228.3	exon20	c.2369C>T	p.T790M	NA	7	55249071	C	T	DP=4;I16=1,0,1,0,69,4761,71,5041,50,2500,42,1764,12,144,25,625;QS=0.543478,0.456522;SGB=-0.379885;RPB=1;MQB=1;BQB=1;MQ0F=0
```


# 2 `src/Somtic_RemoveSite.R`   
## 2.1 依赖:   
`optparse`, `readxl`, `magrittr`, `stringr`, `tidyr`, `Biostrings`包   


## 2.2 用法:   
首次可运行 `Rscript Somtic_RemoveSite.R -h `自动安装依赖的包   
显示如下结果表示依赖安装完成:   

```
  optparse     readxl   magrittr    stringr      tidyr Biostrings 
      TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
Usage: src/Somtic_RemoveSite.R [options]


Options:
	-f FILE, --file=FILE
		

	-s HOTSPOT, --hotspot=HOTSPOT
		

	-d DIR, --dir=DIR
		

	-o OUTDIR, --outdir=OUTDIR
		

	-h, --help
		Show this help message and exit

```
## 2.3 输入/输出文件   
- 输入: **`xxx.hotspot.xls`** , **`xxx.annotate.merge.filter.xls`**
- 输出: **`xxx.annotate.merge.filter.final.xls`** 


## 2.4 输出文件格式   
```
#Sample	Chromosome	Start	End	Ref	Alt	cosmic_id	cosmic_FATHMM_prediction	cosmic_FATHMM_prediction_score	Gene	Region	Classification	CDS_change	Protein_change	Exon_number	Ensembl_transcriptID	Mutant_frequency	Effective_depth	Effective_mutant_depth	Confidence	HotSpot	Oncokb_Mutation_Effect	Oncokb_Oncogenicity	Oncokb_PMIDs_for_Mutation_Effect	CKB_Protein_effect	CKB_Variant_description	1000G_all	Esp6500_all	ExAC_EAS	CytoBand	GenomicSuperDups	CLINSIG	CLNDBN	CLNACC	CLNDSDB	CLNDSDBID	SIFT_prediction	Polyphen2_prediction	MutationAssessor_prediction
Ct1809_1113_NCCL1817_B	1	65312344	65312344	G	T	.	.	.	JAK1	exonic	nonsynonymous SNV	c.C1975A	p.R659S	exon14	ENST00000342505	0.0029	699	2No	No	.	.	.	.	.	.	.	.	1p31.3	.	.	.	.	.	.	T;0.260	P;0.616	L;0.288
Ct1809_1113_NCCL1817_B	1	156843598	156843598	C	T	.	.	.	NTRK1	exonic	nonsynonymous SNV	c.C1024T	p.R342W	exon8	ENST00000524377	0.0024	12503No	No	激活突变	可能致病	1695324, 11313867, 23636398, 22460905	.	.	.	7.7e-05	.	1q23.1	.	.	.	.	.	.	T;0.240	P;0.624	L;0.380
Ct1809_1113_NCCL1817_B	1	156845877	156845877	C	T	.	.	.	NTRK1	exonic	nonsynonymous SNV	c.C1507T	p.H503Y	exon13	ENST00000524377	0.0019	10312No	No	.	.	.	.	.	.	.	.	1q23.1	.	.	.	.	.	.	D;0.721	B;0.425	L;0.450
Ct1809_1113_NCCL1817_B	1	156849086	156849086	C	G	.	.	.	NTRK1	exonic	nonsynonymous SNV	c.C1978G	p.Q660E	exon15	ENST00000524377	0.0017	11542No	No	.	.	.	.	.	.	.	.	1q23.1	.	.	.	.	.	.	T;0.010	B;0.104	N;0.017
Ct1809_1113_NCCL1817_B	3	10183646	10183646	G	A	COSM14370;OCCURENCE=4(kidney)	.	.	VHL	exonic	nonsynonymous SNV	c.G115A	p.G39S	exon1	ENST00000256474	0.0014	1414	2	No	No	.	.	.	.	.	.	.	.	3p25.3	.	.	.	.	.	.	T;0.287	B;0.112	N;0.181
Ct1809_1113_NCCL1817_B	3	10183766	10183766	C	A	.	.	.	VHL	exonic	nonsynonymous SNV	c.C235A	p.R79S	exon1	ENST00000256474	0.0014	1401	2	No	No	.	.	.	.	.	.	.	.	3p25.3	.	.	.	.	.	.	D;0.555	P;0.651	L;0.453
Ct1809_1113_NCCL1817_B	3	10188287	10188287	G	A	COSM1035874;OCCURENCE=1(peritoneum),1(endometrium)	.	.	VHL	exonic	nonsynonymous SNV	c.G430A	p.G144R	exon2	ENST00000256474	0.0016	1229	2	No	No	.	.	.	.	.	.	.	.	3p25.3	.	.	.	.	.	.	.;.	D;0.916	M;0.649
Ct1809_1113_NCCL1817_B	4	1806143	1806143	C	A	.	.	.	FGFR3	exonic	nonsynonymous SNV	c.C1162A	p.L390M	exon9	ENST00000440486	0.0027	752	2	No	No	.	.	.	.	.	.	.	.	4p16.3	.	.	.	.	.	.	T;0.241	B;0.348	L;0.246
```

# 3 RemoveSite_and_db/   

- `Somtic_File_RemoveSite.xlsx`  需要移除的位点
- `Oncokb-annotate.xlsx`   Oncokb注释信息
- `CKB-anV1.xlsx`     CKB注释信息
# ctDNA_output_filter
