# Andres Barboza
# 12/20/2023
# RetrofindR

setwd("~/Documents/GitHub/retrofindR")

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
file.show("test.fasta")

sink("data/query.fasta")
cat(">query\n")
cat(query.seq)
file.show("query.fasta")

