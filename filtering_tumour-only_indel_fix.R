library(tidyverse)
library(RColorBrewer)
library(ggplot2)

#######################################################################################################

# Steps for analysing the data 

# make three folders 1. data 2. output 3. results 

# in the data folder make two more folder of paired and tumour-only data and each folder create vcf and veps
# folder

# run shell script in the vcf folder using vcf files to get the 8 column files 

# chr	pos	ref_allele	alt_allele	ref_reads	var_reads	total_depth	allele_freq


########################################################################################################
getwd()

# for paired data 

#path = "data/paired/vcfs/output/"
#file.names <- dir(path, pattern =".tsv")
#file.names
#mt2_vcf_mpileup <- NULL
#for(i in 1:length(file.names)){
#  b_name <- str_split_fixed(string = basename(file.names[i]), pattern = "[.]",5)[3]
#  print(paste0("* Running analysis for the :" , b_name))
#  vcf_file <- read_tsv(file.path(path, paste0(file.names[i])), col_names = T, col_types = cols(.default=col_character()))
#  vcf_file <- vcf_file %>%
#    mutate_at(vars(ref_reads , var_reads , total_depth ,allel_freq ),funs(as.numeric)) %>%
#    mutate(sample_name = b_name,
#          Location=paste0(chr,":",pos))
#  mt2_vcf_mpileup <- rbind(mt2_vcf_mpileup, vcf_file)
#}

#head(mt2_vcf_mpileup)
#dim(mt2_vcf_mpileup)

# mt2 veps 

#path = "data/paired/veps/0"
#file.names <- dir(path, pattern =".tsv")
#file.names
#mt2_vep_mpileup <- NULL
#for(i in 1:length(file.names)){
#  b_name <- str_split_fixed(string = basename(file.names[i]), pattern = "[.]",5)[3]
#  print(paste0("* Running analysis for the :" , b_name))
#  vep_file <- read_tsv(file.path(path, paste0(file.names[i])), col_names = T, col_types = cols(.default=col_character()),comment = "##")
#  vep_file <- vep_file %>%
#    mutate_at(vars(AFR_AF,ExAC_AF_AFR,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF),funs(as.numeric)) %>%
#    mutate(sample=b_name) %>%
#    select(sample,
#           SYMBOL,
#           Location,
#           Allele,
#          Consequence,
#           VARIANT_CLASS,
#           Feature,
#           STRAND,
#           Protein_position, 
#           Amino_acids, 
#           REF_ALLELE, 
#           BIOTYPE, 
#           SIFT,
#           PolyPhen, 
#           HGVSp, 
#           AFR_AF, 
#           gnomAD_AFR_AF,
#           gnomAD_AMR_AF,
#           gnomAD_ASJ_AF,
#           gnomAD_EAS_AF,
#           gnomAD_FIN_AF,
#           gnomAD_NFE_AF,
#           gnomAD_OTH_AF,
#           gnomAD_SAS_AF,
#           ExAC_AF_AFR)
#  mt2_vep_mpileup <- rbind(mt2_vep_mpileup, vep_file)
#}

#head(mt2_vep_mpileup)
#dim(mt2_vep_mpileup)

# combine vcf and vep files information

#inner_join(mt2_vcf_mpileup, mt2_vep_mpileup, by=c("Location"="Location", 
#                                                  "sample_name"="sample",
#                                                  "ref_allele"="REF_ALLELE",
#                                                  "alt_allele"="Allele")) %>% View()

#mt2_vcf_vep_mpileup <- inner_join(mt2_vcf_mpileup, mt2_vep_mpileup, by=c("Location"="Location",
                                                  # "sample_name"="sample",
                                                  # "ref_allele"="REF_ALLELE",
                                                  # "alt_allele"="Allele"))
#mt2_vcf_vep_mpileup


## for tumour only data 
path = "data/RScriptTroubleshoot/vcfs"
file.names <- dir(path, pattern =".tsv")
file.names
mt2_to_vcf_mpileup <- NULL
for(i in 1:length(file.names)){
  b_name <- str_split_fixed(string = basename(file.names[i]), pattern = "[.]",5)[3]
  print(paste0("* Running analysis for the :" , b_name))
  vcf_file <- read_tsv(file.path(path, paste0(file.names[i])), col_names = T, col_types = cols(.default=col_character()))
  vcf_file <- vcf_file %>%
    mutate_at(vars(ref_reads, var_reads, total_depth, allel_freq),funs(as.numeric)) %>%
    mutate(sample_name = b_name,
           Location=paste0(chr,":",pos))
  mt2_to_vcf_mpileup <- rbind(mt2_to_vcf_mpileup, vcf_file)
}

head(mt2_to_vcf_mpileup)
dim(mt2_to_vcf_mpileup)
#View(mt2_to_vcf_mpileup)

# mt2_to veps
path = "data/RScriptTroubleshoot/veps"
file.names <- dir(path, pattern =".tsv")
file.names
mt2_to_vep_mpileup <- NULL
for(i in 1:length(file.names)){
  b_name <- str_split_fixed(string = basename(file.names[i]), pattern = "[.]",5)[3]
  print(paste0("* Running analysis for the :" , b_name))
  vep_file <- read_tsv(file.path(path, paste0(file.names[i])), col_names = T, col_types = cols(.default=col_character()),comment = "##")
  vep_file <- vep_file %>%
    mutate_at(vars(gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF),funs(as.numeric)) %>%
    mutate(sample=b_name) %>%
    select(sample,
           SYMBOL,
           Location,
           Allele,
           Consequence,
           VARIANT_CLASS,
           Feature,
           STRAND,
           Protein_position, 
           Amino_acids, 
           REF_ALLELE, 
           BIOTYPE, 
           SIFT,
           PolyPhen, 
           HGVSp,
           gnomAD_AFR_AF,
           gnomAD_AMR_AF,
           gnomAD_ASJ_AF,
           gnomAD_EAS_AF,
           gnomAD_FIN_AF,
           gnomAD_NFE_AF,
           gnomAD_OTH_AF,
           gnomAD_SAS_AF)
  mt2_to_vep_mpileup <- rbind(mt2_to_vep_mpileup, vep_file)
}

head(mt2_to_vep_mpileup)
str(mt2_to_vep_mpileup)
dim(mt2_to_vep_mpileup)
table(mt2_to_vep_mpileup$VARIANT_CLASS) ## Still here!


# combine vcf and vep files information *** NOTE: This will only output SNPs, Indels will be lost due to differences in the Location column.
x <- inner_join(mt2_to_vcf_mpileup, mt2_to_vep_mpileup, by=c("Location"="Location", ## The locations vary!
                                                  "sample_name"="sample",
                                                  "ref_allele"="REF_ALLELE",
                                                  "alt_allele"="Allele"))


# Combine SNP VCF and VEP
mt2_to_vcf_vep_mpileup <- inner_join(mt2_to_vcf_mpileup, mt2_to_vep_mpileup, by=c("Location"="Location", 
                                                                         "sample_name"="sample",
                                                                         "ref_allele"="REF_ALLELE",
                                                                         "alt_allele"="Allele"))
View(mt2_to_vcf_vep_mpileup)

# Deal with Indels

# Create VEP and VCFs for missing lines.
anti_mt2_to_vep_mpileup <- anti_join(mt2_to_vep_mpileup, mt2_to_vcf_mpileup, by=c("Location"="Location", 
                                                                                  "sample"="sample_name",
                                                                                  "REF_ALLELE"="ref_allele",
                                                                                  "Allele"="alt_allele"))
anti_mt2_to_vcf_mpileup <- anti_join(mt2_to_vcf_mpileup, mt2_to_vep_mpileup, by=c("Location"="Location", 
                                                                                  "sample_name"="sample",
                                                                                  "ref_allele"="REF_ALLELE",
                                                                                  "alt_allele"="Allele"))

#Check numbers

original_lines <- 227043
combined_lines <- 184275
missing_lines <- original_lines - combined_lines


anti_mt2_to_vep_mpileup[c("LocStart","LocEnd")] <- str_split_fixed(anti_mt2_to_vep_mpileup$Location, "-", 2)

# Rejoin attempt 1 - based on location (StartLoc)

fix_mt2_to_vcf_vep_mpileup <- inner_join(anti_mt2_to_vcf_mpileup, anti_mt2_to_vep_mpileup, by=c("Location"="LocStart", 
                                                                                  "sample_name"="sample",
                                                                                  "ref_allele"="REF_ALLELE",
                                                                                  "alt_allele"="Allele"))


# 

antifix_vcf <- anti_join(anti_mt2_to_vcf_mpileup, anti_mt2_to_vep_mpileup, by=c("Location"="LocStart", 
                                                                                                "sample_name"="sample",
                                                                                                "ref_allele"="REF_ALLELE",
                                                                                                "alt_allele"="Allele"))

antifix_vep <- anti_join(anti_mt2_to_vep_mpileup,anti_mt2_to_vcf_mpileup, by=c("LocStart"="Location", 
                                                                                "sample"="sample_name",
                                                                                "REF_ALLELE"="ref_allele",
                                                                                "Allele"="alt_allele"))

comp <- cbind(antifix_vcf,antifix_vep$LocStart)
comp[c('VCFChr', 'VCFLoc')] <- str_split_fixed(comp$Location, ':', 2)
comp[c('VEPChr', 'VEPLoc')] <- str_split_fixed(comp$`antifix_vep$LocStart`, ':', 2)
comp$VCFLoc <- as.numeric(comp$VCFLoc)
comp$VEPLoc <- as.numeric(comp$VEPLoc)

comp$LocDistance <- comp$VEPLoc - comp$VCFLoc

# This shows that 10473 are in the same location, and the rest are only 1bp off 
# This means if we merge here based on location and sample name only we should be ok, but we'll need to check each time.
table(comp$LocDistance)

# Fix part 2 - combine with location but leave out alleles.
# Minus 1 to VEPLoc in those 7580 that need it.


#Re-built Location column

anti2_vcf <- comp

names(anti2_vcf)[names(anti2_vcf) == 'antifix_vep$LocStart'] <- 'VEPLocationDummy'
anti2_vcf[14307,]
antifix_vep[14307,]
antifix2_mt2_to_vcf_vep_mpileup <- inner_join(anti2_vcf, antifix_vep, by=c("VEPLocationDummy"="LocStart", 
                                                                              "sample_name"="sample"))
### Pay close attention to the warning after the previous line! ###
### It shows where duplicate rows are found and that indel will need manually checking. ###

# In this case the issue is an indel at 14307
# Almost... Successful combination. Re-merge with pileup

dim(antifix2_mt2_to_vcf_vep_mpileup)
colnames(antifix2_mt2_to_vcf_vep_mpileup)
dim(mt2_to_vcf_vep_mpileup)
colnames(mt2_to_vcf_vep_mpileup)
keepnames <- colnames(mt2_to_vcf_vep_mpileup)
dropnames <- c("Location.y","LocEnd","LocDistance","VEPLoc","VCFChr","VCFLoc","VEPChr","VEPLocationDummy","Allele","REF_ALLELE")
antifix2_mt2_to_vcf_vep_mpileup <- antifix2_mt2_to_vcf_vep_mpileup[,!names(antifix2_mt2_to_vcf_vep_mpileup) %in% dropnames]
names(antifix2_mt2_to_vcf_vep_mpileup)[names(antifix2_mt2_to_vcf_vep_mpileup) == "Location.x"] <- "Location"

#There are 4 extra rows in antifix2. 


combined <- rbind(antifix2_mt2_to_vcf_vep_mpileup, mt2_to_vcf_vep_mpileup)
#Now add in the first fix (fix_mt2_to_vcf_vep_mpileup)
#Need to fix columns...
colnames(fix_mt2_to_vcf_vep_mpileup)
colnames(combined)

fix_mt2_to_vcf_vep_mpileup <- fix_mt2_to_vcf_vep_mpileup[,!names(fix_mt2_to_vcf_vep_mpileup) %in% c("Location.y","LocEnd")]

combined2 <- rbind(combined, fix_mt2_to_vcf_vep_mpileup)

mt2_to_vcf_vep_mpileup <- combined2

# remove the tumour only cases that are analysed as mutect2 paired

#mt2_to_vcf_vep_mpileup <- mt2_to_vcf_vep_mpileup[!mt2_to_vcf_vep_mpileup$sample_name %in% unique(mt2_to_vcf_mpileup$sample_name), ]
# sanity check 
#length(unique(mt2_to_vcf_vep_mpileup$sample_name))  # should be 71

# combine paired and tumour only data 

# names(mt2_to_vcf_vep_mpileup)==names(mt2_to_vcf_vep_mpileup) # all true?

#vcf_vep_mpileup <- rbind(mt2_vcf_vep_mpileup, mt2_to_vcf_vep_mpileup)
#vcf_vep_mpileup <-  mt2_to_vcf_vep_mpileup

# vcf_vep_mpileup <- vcf_vep_mpileup %>%
#   filter(sample_name=="B298")

names(mt2_to_vcf_vep_mpileup)
table(mt2_to_vcf_vep_mpileup$VARIANT_CLASS)

# vcf and vep mpileup cleanup e.g. amino acid change and converting consequences from sequencing ontologies to variant classification  
### NOTE! The next command removes secondary values in the Consequences/Variant_Classification column with multiple entries. 
### e.g. "missense_variant,splice_region_variant" will become "Missense_Mutation" and the splice_region_variant info will be lost. 

vcf_vep_mpileup <- mt2_to_vcf_vep_mpileup %>%
  mutate(HGVSp=gsub(".*:","",HGVSp),
         HGVSp=str_replace_all(HGVSp,c("Ala"="A","Arg"="R","Asn"="N","Asp"="D","Asx"="B","Cys"="C","Glu"="E","Gln"="Q","Glx"="Z","Gly"="G","His"="H","Ile"="I","Leu"="L","Lys"="K","Met"="M","Phe"="F","Pro"="P","Ser"="S","Thr"="T","Trp"="W","Tyr"="Y","Val"="V","Xxx"="X","Ter"="*")),
         VARIANT_CLASSIFICATION=gsub("(.*?),.*","\\1",mt2_to_vcf_vep_mpileup$Consequence),
         VARIANT_CLASSIFICATION=if_else(VARIANT_CLASSIFICATION %in% c("missense_variant","coding_sequence_variant","conservative_missense_variant","rare_amino_acid_variant"), "Missense_Mutation",
                                        if_else(VARIANT_CLASSIFICATION%in% c("transcript_amplification","intron_variant","INTRAGENIC","intragenic_variant"),"Intron",
                                                if_else(VARIANT_CLASSIFICATION%in% c("splice_acceptor_variant","splice_donor_variant","transcript_ablation","exon_loss_variant"),"Splice_Site",
                                                        if_else(VARIANT_CLASSIFICATION %in% c("frameshift_variant","protein_altering_variant") & VARIANT_CLASS=="deletion","Frame_Shift_Del",
                                                                if_else(VARIANT_CLASSIFICATION %in% c("frameshift_variant","protein_altering_variant") & VARIANT_CLASS=="insertion","Frame_Shift_Ins",
                                                                        if_else(VARIANT_CLASSIFICATION %in% c("inframe_insertion","disruptive_inframe_insertion","protein_altering_variant") & VARIANT_CLASS=="INS","In_Frame_Ins",
                                                                                if_else(VARIANT_CLASSIFICATION %in% c("inframe_insertion","disruptive_inframe_insertion","protein_altering_variant") & VARIANT_CLASS=="DEL","In_Frame_Del",
                                                                                        if_else(VARIANT_CLASSIFICATION %in% c("incomplete_terminal_codon_variant","synonymous_variant","stop_retained_variant","NMD_transcript_variant"),"Silent",
                                                                                                if_else(VARIANT_CLASSIFICATION %in% c("mature_miRNA_variant","exon_variant","non_coding_exon_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","nc_transcript_variant"),"RNA",
                                                                                                        if_else(VARIANT_CLASSIFICATION %in% c("TF_binding_site_variant","regulatory_region_variant","regulatory_region","intergenic_variant","intergenic_region"), "IGR",
                                                                                                                if_else(VARIANT_CLASSIFICATION %in% c("5_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant"),"5'UTR",
                                                                                                                        if_else(VARIANT_CLASSIFICATION %in% c("initiator_codon_variant","start_lost"),"Translation_Start_Site",
                                                                                                                                if_else(VARIANT_CLASSIFICATION =="stop_gained", "Nonsense_Mutation",
                                                                                                                                        if_else(VARIANT_CLASSIFICATION=="stop_lost","Nonstop_Mutation",
                                                                                                                                                if_else(VARIANT_CLASSIFICATION=="3_prime_UTR_variant", "3'UTR",
                                                                                                                                                        if_else(VARIANT_CLASSIFICATION=="splice_region_variant","Splice_Region",
                                                                                                                                                                if_else(VARIANT_CLASSIFICATION=="upstream_gene_variant","5'Flank","3'Flank"))))))))))))))))),
         GENOME="GRCh37",
         FILTER=if_else(is.na(gnomAD_AFR_AF)| is.na(gnomAD_AMR_AF) | is.na(gnomAD_ASJ_AF)| is.na(gnomAD_EAS_AF)| is.na(gnomAD_FIN_AF) | is.na(gnomAD_NFE_AF)| is.na(gnomAD_OTH_AF) | is.na(gnomAD_SAS_AF), "PASS",
                        if_else (gnomAD_AFR_AF>0.04 | gnomAD_AMR_AF>0.04 | gnomAD_ASJ_AF>0.04 | gnomAD_EAS_AF>0.04 | gnomAD_FIN_AF>0.04 | gnomAD_NFE_AF>0.04 | gnomAD_OTH_AF>0.04 | gnomAD_SAS_AF>0.04, "common_variant","PASS"))) %>%
  select(SAMPLE=sample_name,
         SYMBOL,
         GENOME, 
         CHR=chr,
         START_POS=pos,
         END_POS=pos,
         STRAND,
         VARIANT_CLASSIFICATION,
         VARIANT_TYPE=VARIANT_CLASS,
         REF_ALLELE=ref_allele,
         ALT_ALLELE=alt_allele,
         AA_CHANGE=HGVSp,
         TRANSCRIPT_ID=Feature,
         BIOTYPE,
         SIFT,
         PolyPhen,
         gnomAD_AFR_AF,
         gnomAD_AMR_AF,
         gnomAD_ASJ_AF,
         gnomAD_EAS_AF,
         gnomAD_FIN_AF,
         gnomAD_NFE_AF,
         gnomAD_OTH_AF,
         gnomAD_SAS_AF,
         FILTER,
         TUMOUR_REF_READ=ref_reads,
         TUMOUR_ALT_READ=var_reads,
         TUMOUR_TOTAL_DEPTH=total_depth,
         TUMOUR_VAF=allel_freq)

View(vcf_vep_mpileup)

write_tsv(vcf_vep_mpileup, file.path("Data/RScriptTroubleshoot/results/",paste0("01_VariantsTable_VCF_VEP_RAW.tsv")))


# filtering 

# Constants ---------------------------------------------------------------
vclass <- list()

vclass$genic <- 
  c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
    "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", 
    "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "5'UTR", 
    "3'Flank", "5'Flank", "Intron")

vclass$silent <- c("3'UTR", "5'UTR", "3'Flank", "5'Flank", "Intron", "Silent")

vclass$nonsyn <- setdiff(vclass$genic, vclass$silent)


vcf_vep_mpileup_filtered <- vcf_vep_mpileup %>%
  filter(BIOTYPE=="protein_coding",
         VARIANT_CLASSIFICATION=="Splice_Site" | AA_CHANGE !="-",
         TUMOUR_VAF > 10,
         TUMOUR_ALT_READ >= 4,
         TUMOUR_TOTAL_DEPTH >= 14,
         FILTER=="PASS",
         is.na(gnomAD_NFE_AF) | gnomAD_NFE_AF < 0.01)

head(vcf_vep_mpileup_filtered)
dim(vcf_vep_mpileup_filtered)
write_tsv(vcf_vep_mpileup_filtered, file.path("data/rscripttroubleshoot/results/",paste0("02_VariantsTable_WES_Filtered_By_gnomAD.tsv")))


dim(vcf_vep_mpileup_filtered)

vcf_vep_mpileup_filtered_nonsyn <- vcf_vep_mpileup_filtered[vcf_vep_mpileup_filtered$VARIANT_CLASSIFICATION %in% vclass$nonsyn,]
vcf_vep_mpileup_filtered_silent <- vcf_vep_mpileup_filtered[vcf_vep_mpileup_filtered$VARIANT_CLASSIFICATION %in% vclass$silent,]

## writing output 
dim(vcf_vep_mpileup_filtered_nonsyn)
dim(vcf_vep_mpileup_filtered_silent)
write_tsv(vcf_vep_mpileup_filtered_nonsyn, file.path("data/rscripttroubleshoot/results/",paste0("03_VariantsTable_WES_Filtered_For_Nonsyn_Only.tsv")))
write_tsv(vcf_vep_mpileup_filtered_silent, file.path("data/rscripttroubleshoot/results/", paste0("04_VariantsTable_WES_Filtered_For_Silent_Only.tsv")))

#save(vcf_vep_mpileup, vcf_vep_mpileup_silent, vcf_vep_mpileup_nonsyn, file = "mutations.RData")

FLAGS <- c("MUC4","TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
FLAGS # should display 55 genes
write_tsv(data.frame(FLAGS), file.path("data/rscripttroubleshoot/results/",paste0("55_Flags.tsv")))

#vcf_vep_mpileup_filtered_nonsyn[!vcf_vep_mpileup_filtered_nonsyn$SYMBOL %in% FLAGS,] %>% View()
vcf_vep_mpileup_filtered_nonsyn_Flags_Filtered <- vcf_vep_mpileup_filtered_nonsyn[!vcf_vep_mpileup_filtered_nonsyn$SYMBOL %in% FLAGS,]
dim(vcf_vep_mpileup_filtered_nonsyn_Flags_Filtered)
write_tsv(vcf_vep_mpileup_filtered_nonsyn_Flags_Filtered, file.path("data/rscripttroubleshoot/results/",paste0("05_VariantsTable_WES_Filtered_For_Nonsyn_Flags_Filtered.tsv")))
