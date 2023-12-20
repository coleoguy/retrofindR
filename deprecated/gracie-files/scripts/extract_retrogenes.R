setwd("~/BLAST_Practice/BLAST_Output/")
#load BLAST
blast_output_dris<-read.delim("dris_eejct_fasta_testblast.out", header=F, comment.char = "#")
#assign col names
colnames(blast_output_dris)<-c("eejct","seq","pident","sstart","send","jctstart",
                              "jctend", "align_length","gaps","evalue")
#subset to alignment lengths greater than or equal to #
blast_output_dris<-subset(blast_output_dris, align_length >=35)
eejcts_list<-list()
#loop to remove unique # from geneIDs (FBgn000286_1, _2, _3, etc)
for (i in 1:nrow(blast_output_dris)) {

  eejcts<- strsplit(blast_output_dris$eejct[i], split = "_", fixed=T)[[1]][1]
  
  blast_output_dris$eejct[i] <- eejcts
  
}
#turn blast output into table, all duplicates of gene names are combined and 
#reflected in a frequency #
eejct_pergene<-as.data.frame(table(blast_output_dris$eejct))

eejct.pergene.match <- data.frame(NULL)
#loop through table and match #of blast outputted genes w/sorted alignment 
#score to # of exons in eejct_file for gene 
for(i in 1:nrow(eejct_pergene)){
  for(j in 1:nrow(eejct_list)){
   
    if(eejct_pergene$Var1[i] != eejct_list[[j]]$curgene.name){
      next
  } else {
    if(eejct_pergene$Freq[i] == length(eejct_list[[j]])){
      eejct.pergene.match <- rbind(eejct.pergene.match,
                                   eejct_pergene[i,])
    }
  }
  }
}

#save genes with matched alignments&exons to file on laptop
save(eejct.pergene.match, file="pergenematch")