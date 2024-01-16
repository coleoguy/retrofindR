# Andres Barboza
# 12/20/2023
# RetrofindR

setwd("~/Documents/GitHub/retrofindR")

## FUNCTIONS
# create blast local database
createDB <- function(database.name, subject.fasta) {
  database.command <- paste("makeblastdb -in", subject.fasta, "-out",
                            database.name, "-dbtype 'nucl' -hash_index")
  system(command = database.command)
}
# blastn and load raw results
blastn <- function(database.name, query.file, output.file) {
  blastn.command <- paste("blastn -task blastn -query", query.file, "-db", database.name,
                          "-out data/test_blast_no_header.csv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 10")
  system(blastn.command)
  header <- "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
  header.file <- tempfile(fileext = ".txt")
  writeLines(header, header.file)
  blastn.content <- readLines("data/test_blast_no_header.csv")
  output.content <- c(header, blastn.content)
  writeLines(output.content, output.file)
  blast <- read.csv(output.file, sep = "\t")
  return(blast)
}
# filter blastn results 
filtering <- function(blast, adjacent_homology_cutoff = 40,
                      percent_ident_cutoff = 50,
                      alignment_length_cutoff = 50) {
  pairs_found <- data.frame(
    Value=c(0,0)
  )
  counter <- 0
  for (i in 1:nrow(blast)) {
    for (j in 1:nrow(blast)) {
      if (i != j) {  # Avoid comparing the same row
        diff_value = abs(blast$sstart[i] - blast$send[j])
        if (diff_value <= adjacent_homology_cutoff) {
          counter <- counter + 1
          pairs_found[[paste("value", counter)]] <- c(i, j)
        }
      }
    }
  }
  filtered.blast <- data.frame()
  filtered.blast <- blast[(blast$length >= alignment_length_cutoff),]
  if(length(pairs_found)>1) {
    for (i in 2:length(pairs_found)) {
      if ((blast$pident[as.numeric(pairs_found[[i]][[1]])] + blast$pident[as.numeric(pairs_found[[i]][[2]])])/2 > percent_ident_cutoff){
        filtered.blast <- rbind(filtered.blast, blast[as.numeric(pairs_found[[i]][[1]]),])
        filtered.blast <- rbind(filtered.blast, blast[as.numeric(pairs_found[[i]][[2]]),])
      }
    }
    pairs_found <- pairs_found[,2:length(pairs_found)]
  }
  filtered.blast <- filtered.blast[!duplicated(filtered.blast),]
  filtered.blast <- filtered.blast[(filtered.blast$length - filtered.blast$mismatch - filtered.blast$gapopen) >= alignment_length_cutoff,]
  return(filtered.blast)
}


## ANALYSIS
createDB(database.name = "data/databases/d_melanogaster",
         subject.fasta = "data/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa")


blast <- blastn(database.name = "data/databases/d_melanogaster",
                query.file = "data/fasta/Klp10A.fasta",
                output.file = "data/CSVs/Klp10A_blast.csv")

filtered.blast <- filtering(blast = blast)



##retrocopy screen best practice (from Marques et al.) [https://doi.org/10.1371/journal.pbio.0030357]
#merged adjacent homology matches (distance < 40bp)
#minimum length 50 aminoacids
#amino acid identity < 50%
#aligned over 70% of their length
#merged sequence




### Deprecated #####
# ### Create testing data set 
# exon1 <- sample(c("A","T","C","G"), 120, replace = T)
# exon2 <- sample(c("A","T","C","G"), 75, replace = T)
# exon3 <- sample(c("A","T","C","G"), 48, replace = T)
# 
# raw.seq <- sample(c("A","T","C","G"), 2000, replace = T)
# 
# raw.seq[50:169] <- exon1
# raw.seq[200:274] <- exon2
# raw.seq[400:447] <- exon3
# 
# raw.seq[500:619] <- exon1
# raw.seq[620:694] <- exon2
# raw.seq[695:742] <- exon3
# 
# cds.seq <- raw.seq[500:742]
# 
# test.seq <- paste(unlist(raw.seq), collapse = "")
# query.seq <- paste(unlist(cds.seq), collapse = "")
# 
# sink("data/test.fasta")
# cat(">test\n")
# cat(test.seq)
# #file.show("data/test.fasta")
# sink()
# 
# sink("data/query.fasta")
# cat(">query\n")
# cat(query.seq)
# #file.show("data/query.fasta")
# sink()
#
# 
# 
# #### Test Toy Analysis 
# # create blast local database
# subject.fasta <- "data/test.fasta"
# database.name <- "data/test"
# database.command <- paste("makeblastdb -in", subject.fasta, "-out",
#                           database.name, "-dbtype 'nucl' -hash_index")
# system(command = database.command)
# 
# # do blastn
# query.file <- "data/query.fasta"
# output.file <- "data/test_blast.csv"
# blastn.command <- paste("blastn -query", query.file, "-db", database.name,
#                         "-out data/test_blast_no_header.csv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'")
# system(blastn.command)
# header <- "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
# header.file <- tempfile(fileext = ".txt")
# writeLines(header, header.file)
# blastn.content <- readLines("data/test_blast_no_header.csv")
# output.content <- c(header, blastn.content)
# writeLines(output.content, output.file)
# 
# #read in results
# blast <- read.csv(output.file, sep = "\t")
# #####

