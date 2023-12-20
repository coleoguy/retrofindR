gene_name_true<-read.csv(file="DrosophilaPositionTable.csv", header = TRUE, sep = ",",comment.char = "#")
gene_name_true[gene_name_true$Gene_ID==eejct.pergene.match$Var1,]


#view the number of times a gene appears in the last output by name
view<-blast_output_dris[blast_output_dris$eejct=="FBgn0003346",]
#loop through eejct list and pull out the number of eejcts a gene 
#was predicted to have by name
for(j in 1:length(eejct_list)){
  
  list.names <- c(list.names,eejct_list[[j]]$curgene.name[1])
  
  if(eejct_list[[j]]$curgene.name[1]=="FBgn0003346"){
    k <- nrow(eejct_list[[j]])
    
    l <- j
  } else {
  next
  }
}
#check loop to see if its on the right geneID
nrow(eejct_list[[l]])
#check if geneID exists in blast or eejct list
list.names <- c()
"FBgn0003346" %in% list.names