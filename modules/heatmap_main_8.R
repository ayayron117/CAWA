heatmap.main.8 <- function (gene, avg_APA) {
  
  main_8_tissues <- c("Intestine","Pharynx", "Neuron", "Hypodermis", "Muscle", "Spermatheca", "Vulva_uterus", "Germline")
  
  col_names_geno_tiss <- {}
  
  for (i in 1:length(main_8_tissues)) {
    col_names_geno_tiss[[i]] <- colnames(avg_APA)[grepl(main_8_tissues[i],colnames(avg_APA))]
  }
  
  names(col_names_geno_tiss) <- main_8_tissues
  
  avg_APA_8_tiss <- data.frame(matrix(nrow = 8, ncol = 9))
  geno_names <- sort(gsub("_.*", "", col_names_geno_tiss[[1]]))
  colnames(avg_APA_8_tiss) <- geno_names[c(1,2,3,4,5,7,6,8,9)]
  rownames(avg_APA_8_tiss) <- main_8_tissues
  
  for (i in 1:length(main_8_tissues)) {
    
    avg_apa_vector <- avg_APA[which(rownames(avg_APA) == gene),
                              sort(col_names_geno_tiss[[main_8_tissues[i]]])]
    
    sorted_avg_apa_vec <- avg_apa_vector[,c(1,2,3,4,5,7,6,8,9)]
    
    colnames(sorted_avg_apa_vec) <- gsub("_.*", "", colnames(sorted_avg_apa_vec))
    rownames(sorted_avg_apa_vec) <- main_8_tissues[i]
    
    avg_APA_8_tiss[i, ] <- sorted_avg_apa_vec
    
  }
  
  avg_APA_8_tiss <- avg_APA_8_tiss[,-which(colnames(avg_APA_8_tiss) == "N2D12.D14")]
  
  avg_APA_8_tiss <- avg_APA_8_tiss[,c("N2D1", "N2D6", "LIPL4D1", "LIPL4D6", "DAF2D1", "DAF2D6", "RSKS1D1", "RSKS1D6")]

  return(avg_APA_8_tiss)
  
}










