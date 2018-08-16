current.directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current.directory)

library(data.table)
library(stringr)


all_info <- fread("1000genomes_samples.tsv")
samples <- fread("ALL.chr1.fam")$V1

sum(samples %in% all_info$`Sample name`)
#2504 samples

all_info <- all_info[`Sample name` %in% samples]
all_info[ `Superpopulation code`=="EUR", .(.N), by=`Population name`]

sum(all_info$`Superpopulation code`=="EUR")
# 503 


#--keep accepts a space/tab-delimited text file with family IDs in the first column
#  and within-family IDs in the second column

nation_to_keep <- function(nation){
  ids <- all_info[`Population name`==nation]$`Sample name`
  df <- data.frame("V1"=ids, "V2" = ids)
  write.table(df, file=paste("samples_by_nation/",nation,".keep", sep=""), 
              col.names = F, row.names = F, quote = F)
  }

nations <- unique(all_info$`Population name`)

sapply(nations, function(x) nation_to_keep(x))

ids <- all_info[`Superpopulation name`=="European"]$`Sample name`
df <- data.frame("V1"=ids, "V2" = ids)
write.table(df, file=paste("samples_by_nation/EUR.keep", sep=""), 
            col.names = F, row.names = F, quote = F)


unique(all_info[all_info$`Superpopulation name`=="European"]$`Population name`)