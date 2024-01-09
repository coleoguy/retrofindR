# Andres Barboza
# 12/20/2023
# RetrofindR

setwd("~/Documents/GitHub/retrofindR")

#### Create testing data set ####
exon1 <- sample(c("A","T","C","G"), 120, replace = T)
exon2 <- sample(c("A","T","C","G"), 75, replace = T)
exon3 <- sample(c("A","T","C","G"), 48, replace = T)

raw.seq <- sample(c("A","T","C","G"), 2000, replace = T)

raw.seq[50:169] <- exon1
raw.seq[200:274] <- exon2
raw.seq[400:447] <- exon3

raw.seq[500:619] <- exon1
raw.seq[620:694] <- exon2
raw.seq[695:742] <- exon3

cds.seq <- raw.seq[500:742]

test.seq <- paste(unlist(raw.seq), collapse = "")
query.seq <- paste(unlist(cds.seq), collapse = "")

sink("data/test.fasta")
cat(">test\n")
cat(test.seq)
#file.show("data/test.fasta")
sink()

sink("data/query.fasta")
cat(">query\n")
cat(query.seq)
#file.show("data/query.fasta")
sink()
#####


#### Test Toy Analysis ####
# create blast local database
subject.fasta <- "data/test.fasta"
database.name <- "data/test"
database.command <- paste("makeblastdb -in", subject.fasta, "-out",
                          database.name, "-dbtype 'nucl' -hash_index")
system(command = database.command)

# do blastn
query.file <- "data/query.fasta"
output.file <- "data/test_blast.csv"
blastn.command <- paste("blastn -query", query.file, "-db", database.name,
                        "-out data/test_blast_no_header.csv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'")
system(blastn.command)
header <- "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
header.file <- tempfile(fileext = ".txt")
writeLines(header, header.file)
blastn.content <- readLines("data/test_blast_no_header.csv")
output.content <- c(header, blastn.content)
writeLines(output.content, output.file)

#read in results
blast <- read.csv(output.file, sep = "\t")
#####


#### Test D. melanogaster Analysis ####
# create blast local database
subject.fasta <- "data/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
database.name <- "data/databases/d_melanogaster"
database.command <- paste("makeblastdb -in", subject.fasta, "-out",
                          database.name, "-dbtype 'nucl' -hash_index")
system(command = database.command)

# do blastn
subject.fasta <- "data/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
database.name <- "data/databases/d_melanogaster"
query.file <- "data/fasta/nmda1.fasta"
output.file <- "data/CSVs/nmda1_blast.csv"
blastn.command <- paste("blastn -task blastn -query", query.file, "-db", database.name,
                        "-out data/test_blast_no_header.csv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 10")
system(blastn.command)
header <- "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
header.file <- tempfile(fileext = ".txt")
writeLines(header, header.file)
blastn.content <- readLines("data/test_blast_no_header.csv")
output.content <- c(header, blastn.content)
writeLines(output.content, output.file)

#read in results
blast <- read.csv(output.file, sep = "\t")
#####



#adjacent homology
pairs_within_40 <- data.frame(
  Value=c(0,0)
)
counter <- 0
for (i in 1:nrow(blast)) {
  for (j in 1:nrow(blast)) {
    if (i != j) {  # Avoid comparing the same row
      diff_value = abs(blast$sstart[i] - blast$send[j])
      if (diff_value <= 40) {
        counter <- counter + 1
        pairs_within_40[[paste("value", counter)]] <- c(i, j)
      }
    }
  }
}


filtered.blast <- data.frame()
filtered.blast <- blast[(blast$length >= 50),]
if(length(pairs_within_40)>1) {
  for (i in 2:length(pairs_within_40)) {
    if ((blast$pident[as.numeric(pairs_within_40[[i]][[1]])] + blast$pident[as.numeric(pairs_within_40[[i]][[2]])])/2 > 50){
      filtered.blast <- rbind(filtered.blast, blast[as.numeric(pairs_within_40[[i]][[1]]),])
      filtered.blast <- rbind(filtered.blast, blast[as.numeric(pairs_within_40[[i]][[2]]),])
    }
  }
  pairs_within_40 <- pairs_within_40[,2:length(pairs_within_40)]
}

filtered.blast <- filtered.blast[!duplicated(filtered.blast),]

filtered.blast <- filtered.blast[(filtered.blast$length - filtered.blast$mismatch - filtered.blast$gapopen) >= 50,]




##retrocopy screen best practice (from Marques et al.) [https://doi.org/10.1371/journal.pbio.0030357]
#merged adjacent homology matches (distance < 40bp)
#minimum length 50 aminoacids
#amino acid identity < 50%
#aligned over 70% of their length
#merged sequence

