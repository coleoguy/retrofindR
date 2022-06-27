setwd("~/BLAST_Practice/BLAST_Output/")

blast_output_dris<-read.delim("dris_eejct_fasta_testblast.out", header=F, comment.char = "#")

colnames(blast_output_dris)<-c("eejct","seq","pident","sstart","send","jctstart",
                              "jctend", "align_length","gaps","evalue")

blast_output_dris<-subset(blast_output_dris, align_length >=45)
eejcts_list<-list()
for (i in 1:nrow(blast_output_dris)) {

  eejcts<- strsplit(blast_output_dris$eejct[i], split = "_", fixed=T)[[1]][1]
  
  blast_output_dris$eejct[i] <- eejcts
  
}

eejct_pergene<-as.data.frame(table(blast_output_dris$eejct))

eejct.pergene.match <- data.frame(NULL)

for(i in 1:nrow(eejct_pergene)){
  for(j in 1:length(eejct_list)){
   
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
