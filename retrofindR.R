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
  
  for(i in 1:dim(exons)[1]){
    #Perform initial operations
    #TODO:if [i]=1 do not run
    start <- exons$stop[i]-30
    stop<- exons$stop[i]
    #TODO:if [i]=dim(exons) do not run
    start2<- exons$start[i]
    stop2 <- exons$start[i] + 30
    
    #bind to output dataframe
    eejct <- rbind(eejct,
                   cbind(start,
                         start2,
                         stop,
                         stop2))
    
  }
}


