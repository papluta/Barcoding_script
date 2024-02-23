## HOW TO IDENTIFY BARCODES USING NCBI ONLINE DATA BASE ##
# 1. Go to nBLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
# 2. Upload the fasta file with all sequences or copy your single sequence in the text box
# 3. Optional: in the results table limit your hits to 10 (Sequences producing significant alignments -> show 10) - this will speed up the download but it's not necessary
# 4. in RID (under Job title): 'Download All' as 'Hit Table (csv)'
# 5. Use this script:

install.packages('traits')
library(traits)

# load your hit table:
library(readr)
library(dplyr)

seq.raw <- read.csv('~/XE40ZAJG013-Alignment-HitTable.csv', header = F)

## FOR ONE SINGLE BEST HIT
best_taxon_fun <- function(x) {
  colnames(x) <- c('query', 'subject', '%identity', 'alignment_length', 'mismatches', 'gap_opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score')
  seq <- distinct(x, `query`, .keep_all = T)
  acc.no <- seq[,2]
  acc.no <- unlist(acc.no)
  taxons <- ncbi_byid(ids=acc.no)[,1]
  seq$taxons <- taxons
  seq$id <- sub(".*_","", seq$`query`)
  return(seq)
  }

## FOR 3 BEST HITS PER SAMPLE
best3_taxons_fun <- function(x) {
  colnames(x) <- c('query', 'acc.n', '%identity', 'alignment_length', 'mismatches', 'gap_opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score')
  seq <- x %>%
    group_by(`query`) %>%
    filter(row_number() <= 3)
  acc.no <- seq[,2]
  acc.no <- unlist(acc.no)
  taxons <- ncbi_byid(ids=acc.no)[,1]
  seq$taxons <- taxons
  seq$id <- sub(".*_","", seq$`query`)
  return(seq)
}

## get the taxons:
results.best <- best_taxon_fun(seq.raw) ## if there's an error or timeout, split your seq.raw into smaller parts of ~200 samples:
#results.best <- lapply(list(seq.raw[1:200,], seq.raw[201:400,], seq.raw[401:600,], seq.raw[601:nrow(seq.raw),]), best_taxon_fun) %>% bind_rows() #each element in the list is raw.seq split by rows

results.3best <- best3_taxons_fun(seq.raw)

# add sample names - if you have a file with the well position on the plate and ID of the sample:
sample <- read.csv('~/Sample_setup.csv') # name the well position in your file 'Plate' so the below code works:

results <- left_join(sample, results.best, by = join_by('Plate' == 'id'))


write.csv(results, file = '~/Barcoding1_results.csv', row.names = F)
