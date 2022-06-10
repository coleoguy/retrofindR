setwd("~/R stuff")
library(seqinr)
genome<-read.fasta("Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa")
load("drosjctlist")


eejct_code_list <- vector("list",length(eejct_list))

for(i in 1:length(eejct_list)){
  cur.table <- eejct_list[[i]]
  curr.seq <- genome[[ which(names(genome) == cur.table$cur.seq[1])]]
  

    for(j in 1:nrow(cur.table)){
      eejct_code_list[[i]][[j]] <- curr.seq[c(cur.table[j,2]:cur.table[j,3],cur.table[j,4]:cur.table[j,5])] 
      
     # eejct_code_list[[length(eejct_code_list) + 1]] <- eejct_code
    }
  
  
}













q_unique <- paste0("seq_",cur.table$cur.seq[1])

seq_code_list[[length(seq_code_list) + 1]] <- assign(seq_unique,eejct_code)




if(eejct_list[[i]]$cur.seq[j] == names(genome)[j]){
  curr_seq <- genome[[j]]
}






eejct_seqs_list[[length(master_eejct_code_list) + 1]] <- eejct_code_list