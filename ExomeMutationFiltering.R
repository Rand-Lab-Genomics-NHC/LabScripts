library(tidyverse)
library(ggplot2)
library(maftools)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
# MT2 
# vcf

setwd("D:\\Bioinfo_Analysis/dlbcl_cohort/ExomeMutationFiltering/MT2/vcf/")

# reads .tsv files format: chr    pos   alt  ref_reads       var_reads   allele_freq


path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2 <- NULL
for(i in 1:length(file.names)){
  data <- read.table(file.names[i], sep = "\t", header = FALSE,as.is = TRUE)
  colnames(data) <- c("chr", "pos","alt" ,"ref_read", "var_reads", "allele_freq")
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[3]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(data))
  data <- cbind(Sample_Name, data)
  data$Location <- paste0(data$chr,":",data$pos)
  data <- data[,c("Sample_Name", "Location", "chr", "pos","alt", "ref_read", "var_reads", "allele_freq")]
  MT2 <- rbind(MT2, data)
}

dim(MT2)
head(MT2)
table(MT2$Sample_Name)

# MT2 
# vep

setwd("D:\\Bioinfo_Analysis/dlbcl_cohort/ExomeMutationFiltering/MT2/vep/")

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_veps <- NULL
for (i in 1:length(file.names)){
  data <- read_tsv(file.names[i], col_names = TRUE, col_types = cols(gnomAD_NFE_AF= col_double(),Location=col_character()))
  filtered.data <- data %>%
    select("Location", "Allele", "Consequence",	"Protein_position",	"Amino_acids","Codons","STRAND",	"VARIANT_CLASS",	"SYMBOL","BIOTYPE","SIFT","PolyPhen", "gnomAD_NFE_AF") %>%
    mutate(gnomAD_NFE_AF= replace_na(gnomAD_NFE_AF,0)) %>%
    filter(BIOTYPE=="protein_coding", Amino_acids !="-", str_detect(Consequence, "synonymous",negate = TRUE), str_detect(Consequence, "splice_region_variant", negate = TRUE), gnomAD_NFE_AF < 0.01)
  Amino_acid <- str_split_fixed(filtered.data$Amino_acids,pattern = "/",2)
  Protein_position <- str_split_fixed(filtered.data$Protein_position,"/",2)
  AA <- paste0(Amino_acid[,1],Protein_position[,1],Amino_acid[,2])
  filtered.data <- mutate(filtered.data,Amino_acids=AA)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[3]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(filtered.data))
  filtered.data <- cbind(Sample_Name, filtered.data)
  MT2_veps <- rbind(MT2_veps, filtered.data)
}

head(MT2_veps)
dim(MT2_veps)

####################################################################################################################
# tumour_only 
# vcf 
####################################################################################################################
setwd("D:\\Bioinfo_Analysis/dlbcl_cohort/ExomeMutationFiltering/MT2_tumour_only/vcf/")

# reads .tsv files format: chr    pos   alt  ref_reads       var_reads   allele_freq

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_tumour_only <- NULL
for(i in 1:length(file.names)){
  data <- read.table(file.names[i], sep = "\t", header = FALSE, as.is = TRUE)
  colnames(data) <- c("chr", "pos", "alt","ref_read", "var_reads", "allele_freq")
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[1]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(data))
  data <- cbind(Sample_Name, data)
  data$Location <- paste0(data$chr,":",data$pos)
  data <- data[,c("Sample_Name", "Location", "chr", "pos","alt", "ref_read", "var_reads", "allele_freq")]
  MT2_tumour_only <- rbind(MT2_tumour_only, data)
}

dim(MT2_tumour_only)
head(MT2_tumour_only)

setwd("D:\\Bioinfo_Analysis/dlbcl_cohort/ExomeMutationFiltering/MT2_tumour_only/vep/")

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_tumour_only_veps <- NULL
for (i in 1:length(file.names)){
  data <- read_tsv(file.names[i], col_types = cols(gnomAD_NFE_AF= col_double(), Location=col_character()))
  colnames(data)
  filtered.data <- data %>%
    select("Location", "Allele", "Consequence",	"Protein_position",	"Amino_acids","Codons","STRAND",	"VARIANT_CLASS",	"SYMBOL","BIOTYPE","SIFT","PolyPhen", "gnomAD_NFE_AF") %>%
    mutate(gnomAD_NFE_AF= replace_na(gnomAD_NFE_AF,0)) %>%
    filter(BIOTYPE=="protein_coding", Amino_acids !="-", str_detect(Consequence, "synonymous",negate = TRUE),str_detect(Consequence, "splice_region_variant", negate = TRUE), str_detect(VARIANT_CLASS, "substitution", negate = TRUE), gnomAD_NFE_AF < 0.01)
  Amino_acid <- str_split_fixed(filtered.data$Amino_acids,pattern = "/",2)
  Protein_position <- str_split_fixed(filtered.data$Protein_position,"/",2)
  AA <- paste0(Amino_acid[,1],Protein_position[,1],Amino_acid[,2])
  filtered.data <- mutate(filtered.data,Amino_acids=AA)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[1]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(filtered.data))
  filtered.data <- cbind(Sample_Name, filtered.data)
  MT2_tumour_only_veps <- rbind(MT2_tumour_only_veps, filtered.data)
}

head(MT2_tumour_only_veps)
dim(MT2_tumour_only_veps)
table(MT2_tumour_only_veps$Sample_Name)

####################################################################################################################
# vcf and vep joining 
####################################################################################################################

# MT2 files

MT2 #vcf 
MT2_veps #veps

dim(MT2)
dim(MT2_veps)
head(MT2)
head(MT2_veps)

MT2_veps$Location <- gsub("\\-.*","",MT2_veps$Location) # trim after - to include the inframe or framshit variants

MuTect2 <- inner_join(MT2_veps, MT2, by = c("Sample_Name"="Sample_Name","Location"="Location"))
dim(MuTect2)
head(MuTect2)

# MT2_tumour_only files

MT2_tumour_only #vcf 
MT2_tumour_only_veps #veps

dim(MT2_tumour_only)
dim(MT2_tumour_only_veps)
head(MT2_tumour_only)
head(MT2_tumour_only_veps)

MT2_tumour_only_veps$Location <- gsub("\\-.*","",MT2_tumour_only_veps$Location) # trim after - to include the inframe or framshit variants
MuTect2_Tumour_Only <- inner_join(MT2_tumour_only_veps, MT2_tumour_only, by = c("Sample_Name"="Sample_Name","Location"="Location"))
dim(MuTect2_Tumour_Only)
head(MuTect2_Tumour_Only)

# check 
table(colnames(MuTect2)==colnames(MuTect2_Tumour_Only)) # true for correct order

dlbcl_mutations <- rbind(MuTect2,MuTect2_Tumour_Only)
table(order(dlbcl_mutations$SYMBOL))

# Further processing 

dlbcl_mutations <- dlbcl_mutations %>% 
  mutate(Codons=str_replace_all(Codons,"[:lower:]","")) # remove lower letter from the Codon column 
dlbcl_mutations <- dlbcl_mutations %>%
  separate(Codons,c("Ref_allele","Alt_allele"),"/")

dlbcl_mutations <- dlbcl_mutations[dlbcl_mutations$var_reads >= 4 & dlbcl_mutations$ref_read >= 10,] # 
dim(dlbcl_mutations)
head(dlbcl_mutations)


setwd("D:\\Bioinformatics_Analysis/dlbcl_cohort/ExomeMutationFiltering")
write_tsv(dlbcl_mutations, "dlbcl_mutations_GATK4.tsv")

tbl <- as.data.frame(table(MuTect2$SYMBOL))
tbl[order(tbl$Freq,decreasing = TRUE),]



########################################################################
#### plots 
########################################################################
p1<- dlbcl_mutations %>%
  select(Sample_Name) %>%
  group_by(Sample_Name) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(Sample_Name, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Sample Names") +
  ylab("Number of Mutations") +
  ggtitle("Mutations per Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p2<- dlbcl_mutations %>%
  select(chr) %>%
  group_by(chr) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(chr, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Sample Names") +
  ylab("Number of Mutations") +
  ggtitle("Mutations Load by Chromosome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p3<- dlbcl_mutations %>%
  select(SYMBOL) %>%
  group_by(SYMBOL) %>%
  count() %>%
  filter(n>20) %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(SYMBOL, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Genes Names") +
  ylab("Number of Mutations > 20") +
  ggtitle("Mutations Load by Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p4<- dlbcl_mutations %>%
  select(Consequence) %>%
  group_by(Consequence) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(Consequence, n),y=n, fill=Consequence)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  coord_flip() +
  ylab("Variant Counts") +
  ggtitle("Variant Classification") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),axis.title.y = element_blank())

p5<- dlbcl_mutations %>%
  select(VARIANT_CLASS) %>%
  group_by(VARIANT_CLASS) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(VARIANT_CLASS, n),y=n, fill=VARIANT_CLASS)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  coord_flip() +
  ylab("Variant Counts") +
  ggtitle("Variant Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),axis.title.y = element_blank())


mycolors <- colorRampPalette(brewer.pal(8,"Dark2"))(31)

p6<- dlbcl_mutations %>%
  select(Sample_Name, allele_freq) %>%
  group_by(Sample_Name) %>%
  ggplot(aes(x=allele_freq, fill=Sample_Name)) +
  geom_density(alpha= 0.5,position = "stack",show.legend = FALSE)+
  scale_fill_manual(values = mycolors) +
  #facet_wrap(~Sample_Name) +
  xlab("Allelic Frequency") +
  ylab("Density") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),)

png("MutationsProfiles_dlbcl.png",units = "in",width = 18,height = 12,res = 300)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
graphics.off()

dlbcl_mutations %>%
  select(Sample_Name,ref_read, var_reads) %>%
  mutate(Allelic_Depth= ref_read+var_reads) %>%
  ggplot(aes(x=Sample_Name, y = Allelic_Depth, fill=Sample_Name)) +
  geom_boxplot(show.legend = FALSE) +
  #facet_wrap(~Sample_Name) +
  xlab("Sample Names") +
  ylab("Allelic Depth") +
  ggtitle("Allelic Depth per Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

test<- dlbcl_mutations %>%
  group_by(SYMBOL, Sample_Name) %>%
  slice(1) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n) *100) %>%
  complete(SYMBOL, Sample_Name, fill=list(percent =0)) %>%
  select(SYMBOL,Sample_Name, percent)

test<- dlbcl_mutations %>%
  select(SYMBOL, Sample_Name) %>%
  group_by(SYMBOL, Sample_Name) %>%
  slice(1) %>%
  summarise(count=n()) %>%
  spread(Sample_Name, count, fill=0) %>%
  column_to_rownames(var = "SYMBOL") 


test <- rowSums(test[,sapply(test, is.numeric)])  
test <- data.frame(test)
test<- test %>%
  rownames_to_column("SYMBOL") %>%
  mutate(percent = test/31*100)
View(test)

#################### oncoplot
test <- dlbcl_mutations %>%
  select(SYMBOL,Sample_Name,VARIANT_CLASS) %>%
  pivot_wider(names_from = VARIANT_CLASS, values_from = VARIANT_CLASS,values_fill = list(VARIANT_CLASS=0),values_fn = list(VARIANT_CLASS = length))
  

snv <- test %>%
  select(SYMBOL,Sample_Name,SNV) %>%
  pivot_wider(names_from = Sample_Name, values_from = SNV, values_fill = list(SNV=0),values_fn = list(SNV=length)) %>%
  column_to_rownames("SYMBOL")%>%
  as.matrix()

insertion <- test %>%
  select(SYMBOL,Sample_Name,insertion) %>%
  pivot_wider(names_from = Sample_Name, values_from = insertion, values_fill = list(insertion=0),values_fn = list(insertion=length)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
view(insertion)

mat_list <- list(snv,insertion)
names(mat_list) <- c("snv", "insertion")
mat_list 

col = c(snv = "red", insertion = "blue")



col = c("insertion" = "blue", "snv" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["insertion"], col = NA))
  },
   # small green
  snv = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["snv"], col = NA))
  }
)

heatmap_legend_param = list(title = "Alternations", at = c("insertion", "snv"), 
                            labels = c("Insertion", "SNV"))

oncoPrint(mat_list, 
          alter_fun = alter_fun, 
          col = col, heatmap_legend_param = heatmap_legend_param)
