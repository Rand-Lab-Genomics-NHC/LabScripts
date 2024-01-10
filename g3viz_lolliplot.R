# G3Viz Lollipop plot script - Given that GenVisR is deprecated this 
# script is to replace the lolliplot script that relied on GenVisR but using G3Viz instead.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library("BiocManager")
BiocManager::install("g3viz", force=T)
BiocManager::install("GenomicAlignments")
library("g3viz")
setwd("~/Data/Lolliplot_Workings_10012024/")

#maf.file <- read.csv2("Filtered_Nonsyn_Flags.csv", sep=',', header=T)
#maf.file <- system.file("extdata","~/Data/Lolliplot_Workings_10012024/Filtered_Nonsyn_Flags.csv", package='g3viz')

# prefile <- read.csv2("Filtered_Nonsyn_Flags.csv", sep=',', header=T)
# prefile$AA_CHANGE <- gsub("p.", "", prefile$AA_CHANGE)
# write.table(prefile, "Filtered_Nonsyn_Flags_treat.csv", quote=F, sep=',',row.names = F)

maf.file <- "~/Data/Lolliplot_Workings_10012024/Filtered_Nonsyn_Flags.csv"

data <- g3viz::readMAF(maf.file,
                       gene.symbol.col = "SYMBOL",
                       variant.class.col = "VARIANT_CLASSIFICATION",
                       protein.change.col = "AA_CHANGE",
                       if.parse.aa.pos = T,
                       if.parse.mutation.class = T,
                       mutation.type.to.class.df = NA,
                       sep = ",")


chart.options <- g3Lollipop.theme(theme.name = "default",
                                  title.text = "Whatever gene")

g3Lollipop(data,
           gene.symbol = "KIF17",
           plot.options = chart.options,
           protein.change.col = "AA_CHANGE",
           gene.symbol.col = "SYMBOL",
           output.filename = "default_theme")



