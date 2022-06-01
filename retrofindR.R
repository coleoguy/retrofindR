# This is a script to work on the initial steps of the retrofindR package.

gff <- read.delim("Drosophila_melanogaster.BDGP6.32.106.gff3", header = F,
                  comment.char = "#")
colnames(gff) <- c("seq","source","feature","start","stop","score","strand","phase","info")


# we need an exon exon junction file
# here lets create an empty vector that we will fill with the data we need
eejct <- data.frame()

# this loop works its way through the gff file
for(i in 1:nrow(gff)){
  # we need to identify when we are on a row that starts a gene
  if(gff$feature[i] == "gene"){
     x <- strsplit(gff$info[i], split = ";", fixed=T)[[1]][1]
     curgene.name <- strsplit(x, split=":", fixed=T)[[1]][2]
     cur.start<-gff$start[i]
     cur.end<-gff$stop[i]
     cur.seq<-gff$seq[i]
#how to sort by location and sequence
  exons<-subset(gff,gff$start>=cur.start & gff$stop<=cur.end & gff$feature=="exon"& gff$seq==cur.seq)[,c(1,4,5)]
}
exons$gene<-curgene.name
exons<-exons[!duplicated(exons$start),]
}
  # add another if statement that is true when you have two exons 
  # that are in the same gene (curgene) when it is true
  # store the start and stop positions for the part 
  # you want to become your exon exon jct sequences
#if statement printing junctions:
if(nrow("exons")>1){
  i <- 2
  working <- T
  while(working){
    start <- exons$stop[i-1]-30
    stop<- exons$stop[i-1]
    start2<- exons$start[i]
    stop2 <- exons$start[i] + 30
    if(stop2 < stop){
      working <- F
    }else{
      #bind to output dataframe
      eejct <- rbind(eejct,
                     cbind(cur.seq,start,stop,start2,
                           stop2,curgene.name))
    }
    if(i == nrow(exons)){
      working <- F
    }
    i <- i + 1
  }
}

#TODO: run code through gff file and save each eejct data set 
#as one component of a list