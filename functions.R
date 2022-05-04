dotplots <- function(goi){
  
  goi <- read.table(input$user_file$datapath)
  rownames(goi) <- goi[,1]
  
  if(!any(rownames(goi) %in% rownames(bulkseq_mat))){
    return()
  }
  
  
  goi <- goi[which(rownames(goi) %in% rownames(bulkseq_mat)==TRUE),]
  
  mat <- bulkseq_mat[,1:154]
  
  tcounts <- t(mat[goi, ]) %>%
    base::merge(bulkseq_colData, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):(ncol(.)))

  tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol
  
  tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, expression) %>% dplyr::summarise(expression=median(expression))
  tcounts_med[!tcounts_med$Condition %in% "Ipsi", ] #remove all injured samples
  
  g <- ggplot(tcounts_med, aes(x=Population, y = symbol)) #group by population
  g <- g + scale_colour_viridis_c(option = "magma", end = .90)
  g <- g + geom_point(aes(col=expression, size=expression))
  g <- g + theme_bw() + theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size=10, angle = 45, hjust= 1),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank())
  g <- g + scale_x_discrete(labels=c(
    "TBAC" = "A\u03b4-LTMR + A\u03b2-RA-LTMR",
    "CRTH" = "C-LTMR",
    "MRTD" = "NP",
    "CGRT" = "PEP",
    "TDNV" = "Nociceptors"))

  return(g)
} 
