# This is a script to work on the initial steps of the retrofindR package.

gff <- read.delim("Drosophila_melanogaster.BDGP6.32.106.gff3", header = F,
                  comment.char = "#")
colnames(gff) <- c("seq","source","feature","start","stop","score","strand","phase","info")


# we need an exon exon junction file
# here lets create an empty vector that we will fill with the data we need
eejct <- vector()

# this loop works its way through the gff file
for(i in 1:nrow(gff)){
  # we need to identify when we are on a row that starts a gene
  if(gff$feature[i] == "gene"){
     x <- strsplit(gff$info[i], split = ";", fixed=T)[[1]][1]
     curgene.name <- strsplit(x, split=":", fixed=T)[[1]][2]
     cur.start<-gff$start[i]
     cur.end<-gff$stop[i]
     cur.seq<-gff$seq[i]
  }
  
  # add another if statement that is true when you have two exons 
  # that are in the same gene (curgene) when it is true
  # store the start and stop positions for the part 
  # you want to become your exon exon jct sequences
}
#sorting by gene name doesnt work
library(stringr)
gff[str_detect(gff$info,curgene.name)]

#how to sort by location and sequence
exons<-subset(gff,gff$start>=cur.start & gff$stop<=cur.end & gff$feature=="exon"& gff$seq==cur.seq)
