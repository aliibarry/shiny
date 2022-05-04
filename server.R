##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   updated 09/2021                    ##
##   for non-commercial use only        ##
##########################################

#server doc

shinyServer(function(input, output, session) {
  
  load("drg-directory.RData")
  
    ### linking to other tabs
    observeEvent(input$link_to_tabsummary, {
        newvalue <- "tabsummary"
        updateTabItems(session, "tabs", newvalue)
    })
    
    #output$genesearch <- eventReactive(input$goButton,{input$genequery})

    data <- reactive({ 
      req(input$user_file) 
      
      goi <- read.table(input$user_file$datapath)
      rownames(goi) <- goi[,1]
      
      goi <- goi[which(rownames(goi) %in% rownames(bulkseq_mat)==TRUE),]
      
      return(goi)
    })
    
    
    output$goi_table <- DT::renderDataTable({
      
      req(input$user_file);

        datatable <- bulkseq_mat[data(),]

        return(datatable)

    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("shinydata_usergrep", ".csv", sep = "")
      },
      
      content = function(file) {
        
        goi <- data()
        datatable <- bulkseq_mat[goi,]
        
        write.csv(datatable, file, row.names = FALSE)
      })
    
    
    output$bulkseq_dots <- renderPlot({
      
      req(input$user_file);
      
      mat <- bulkseq_mat[,1:154]
      goi <- data()
      
      tcounts <- t(mat[goi, ]) %>%
        base::merge(bulkseq_colData, ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-length(goi)+1):(ncol(.)))
      
      ## add gene symbols for goi terms
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
    })
    
    output$bulkseq_lines <- renderPlot({
        
        req(input$user_file);
      
        mat <- bulkseq_mat[,1:154]
        goi <- data()
        
        tcounts <- t(mat[goi, ]) %>%
            base::merge(bulkseq_colData, ., by="row.names") %>%
            gather(gene, expression, (ncol(.)-length(goi)+1):(ncol(.)))
        
        ## add gene symbols for goi terms
        tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol
        
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, Timepoint, symbol) %>% dplyr::summarise(expression=median(expression))
        
        g <- ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint)) ) 
        g <- g + scale_colour_viridis(discrete=TRUE, end = .80)
        g <- g + geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol))
        g <- g + facet_wrap(~symbol, ncol=3)
        g <- g + theme_bw() + theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(size=10, angle = 45, hjust= 1),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank())

        return(g)
    })
    
    
    
    ### leaflet map for contact details
    output$myMap <- renderLeaflet({
        m <- leaflet() %>% addTiles()
        m <- m %>% setView( -1.2217, 51.76529, zoom = 16)
        m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
    })
    
    #volcanoServer("volcano_module")
    
    
})
