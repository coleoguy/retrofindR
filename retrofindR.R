# This is a script to work on the initial steps of the retrofindR package.

gff <- read.delim("Drosophila_melanogaster.BDGP6.32.106.gff3", header = F,
                  comment.char = "#")
colnames(gff) <- c("seq","source","feature","start","stop","score","strand","phase","info")


# we need an exon exon junction file
# here lets create an empty vector that we will fill with the data we need
eejct <- data.frame()

eejct_list <- list()
# this loop works its way through the gff file
for(i in 1:nrow(gff)){
  eejct <- data.frame()
  #exons <- data.frame()
  
  # we need to identify when we are on a row that starts a gene
  if(gff$feature[i] == "gene"){
     x <- strsplit(gff$info[i], split = ";", fixed=T)[[1]][1]
     curgene.name <- strsplit(x, split=":", fixed=T)[[1]][2]
     cur.start<-gff$start[i]
     cur.end<-gff$stop[i]
     cur.seq<-gff$seq[i]
#how to sort by location and sequence
  exons<-subset(gff,gff$start>=cur.start & gff$stop<=cur.end & gff$feature=="exon"& gff$seq==cur.seq)[,c(1,4,5)]
  
  #Remove dulications and add gene_id
  exons$gene<-curgene.name
  exons<-exons[!duplicated(exons$start),]
  
  if(nrow(exons) > 1){
    j <- 2
    working <- T
    while(working){
      start <- exons$stop[j-1]-30
      stop<- exons$stop[j-1]
      start2<- exons$start[j]
      stop2 <- exons$start[j] + 30
      if(stop2 < stop){
        working <- F
      }else{
        #bind to output dataframe
        eejct <- rbind(eejct,
                       cbind(cur.seq,start,stop,start2,
                             stop2,curgene.name))
      }
      if(j == nrow(exons)){
        working <- F
      }
      j <- j + 1
    }
    
    
    eejct_unique <- paste0("eejct_",eejct$curgene.name[1])
    
    eejct_list[[length(eejct_list) + 1]] <- assign(eejct_unique,eejct)
  }
}
# exons$gene<-curgene.name
# exons<-exons[!duplicated(exons$start),]
# 
#   j <- 2
#   working <- T
#   while(working){
#     start <- exons$stop[j-1]-30
#     stop<- exons$stop[j-1]
#     start2<- exons$start[j]
#     stop2 <- exons$start[j] + 30
#     if(stop2 < stop){
#       working <- F
#     }else{
#       #bind to output dataframe
#       eejct <- rbind(eejct,
#                      cbind(cur.seq,start,stop,start2,
#                            stop2,curgene.name))
#     }
#     if(k == nrow(exons)){
#       working <- F
#     }
#     k <- k + 1
#   }
# 
# 
# eejct_unique <- paste0("eejct_",eejct$curgene.name[1])
# 
# eejct_list[[length(eejct_list) + 1]] <- assign(eejct_unique,eejct)
}

eejct_list

gff$info[9]
  # add another if statement that is true when you have two exons 
  # that are in the same gene (curgene) when it is true
  # store the start and stop positions for the part 
  # you want to become your exon exon jct sequences
#if statement printing junctions:
# if(nrow("exons")>1){
#   i <- 2
#   working <- T
#   while(working){
#     start <- exons$stop[i-1]-30
#     stop<- exons$stop[i-1]
#     start2<- exons$start[i]
#     stop2 <- exons$start[i] + 30
#     if(stop2 < stop){
#       working <- F
#     }else{
#       #bind to output dataframe
#       eejct <- rbind(eejct,
#                      cbind(cur.seq,start,stop,start2,
#                            stop2,curgene.name))
#     }
#     if(i == nrow(exons)){
#       working <- F
#     }
#     i <- i + 1
#   }
# }

#TODO: run code through gff file and save each eejct data set 
#as one component of a list

length(eejct_list)