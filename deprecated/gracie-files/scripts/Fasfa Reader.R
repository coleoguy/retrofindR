setwd("~/R stuff/retrofindR")
#run package
library(seqinr)
#read in FASTA file
genome<-read.fasta("Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa")
#read in eejct list
load("drosjctlist")


eejct_code_list <- list()

gene_names <- c()

#loop through eejct list and match gene names from eejct list to 
#gene names in organism genome FASTA
for(i in 1:length(eejct_list)){
  cur.table <- eejct_list[[i]]
  curr.seq <- genome[[ which(names(genome) == cur.table$cur.seq[1])]]
  
  
    for(j in 1:nrow(cur.table)){
      cur.table$curgene.name[j] <- paste0(cur.table$curgene.name[j],"_",j)
      #paste out sequences for each eejct based on FASTA sequence of 
      #organism and locations from retrofindR.R
      eejct_code_list[[length(eejct_code_list) + 1]] <- curr.seq[c(cur.table[j,2]:cur.table[j,3],cur.table[j,4]:cur.table[j,5])] 
      gene_names <- c(gene_names,cur.table$curgene.name[j])
     #<- eejct_code 
    }
  
  
}

eejct_code_list[[1]]

#Check gene names 
gene_names

#Add names to eejct_code_ist
names(eejct_code_list)<-gene_names

write.fasta(eejct_code_list,names(eejct_code_list), open="w",as.string = FALSE,
            nbchar = 60, file.out = "writefastatest.fasta")






