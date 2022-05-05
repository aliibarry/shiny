##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

#server doc

shinyServer(function(input, output, session) {
  
  load("drg-directory.RData")
  
  updateSelectizeInput(session, 
                       inputId = "geneid", 
                       label = "",
                       choices = bulkseq_mat[,155], 
                       server = TRUE,
                       selected = "Atf3"
                       #options = list(placeholder = 'select a gene', 'plugins' = list('remove_button'))
                       )
  
    ### linking to other tabs
    observeEvent(input$link_to_tabsummary, {
        newvalue <- "tabsummary"
        updateTabItems(session, "tabs", newvalue)
    })
    
    observeEvent(input$link_to_tables, {
      newvalue <- "tabdata"
      updateTabItems(session, "tabs", newvalue)
    })
    
    observeEvent(input$link_to_wald, {
      newvalue <- "tabwald"
      updateTabItems(session, "tabs", newvalue)
    })
    
    #output$genesearch <- eventReactive(input$goButton,{input$genequery})

    # data <- reactive({
    #     req(input$user_file)
    # 
    #     goi <- read.table(input$user_file$datapath)
    #     rownames(goi) <- goi[,1]
    # 
    #     goi <- goi[which(rownames(goi) %in% rownames(bulkseq_mat)==TRUE),]
    # 
    #     return(goi)
    # })
        
    data <- reactive({
      req(input$geneid)

      filteredmat <- bulkseq_mat[bulkseq_mat[,155] %in% input$geneid,]
      #filteredmat <- bulkseq_mat[bulkseq_mat[,155] %in% input$user_file,]
        
      return(filteredmat)
          
      })
    
    usercontrast <- reactive({
      req(input$contrast)
      res <- read.csv(input$contrast)
      return(res)
    })
    
    output$contrast_table <- DT::renderDataTable({
      
      req(usercontrast())
      
      DT::datatable(
        usercontrast(),
        width = 12,
        #style="default",
        #fillContainer = TRUE,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
    
    
    
    
    output$goi_table <- DT::renderDataTable({
        
      req(data())
      
      datatable <- data()
      
      DT::datatable(
        datatable,
        style="default",
        #fillContainer = TRUE,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 5)
      )
          
       # if (input$user_file == TRUE){
       #   datatable <- bulkseq_mat[data(),]
       # }
       #  else {
       #    datatable <- bulkseq_mat[grep(input$geneid, bulkseq_mat[,155]),]
       #  }
       # 
       #  return(datatable)

    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("drgdirectory_search", ".csv", sep = "")
      },
      
      content = function(file) {
        
        datatable <- data()
        write.csv(datatable, file, row.names = FALSE)
      })
    
    PlotHeight = reactive(
      return(length(data()))
    )
    
    output$bulkseq_dots <- renderPlotly({
      
      req(data());
      
      mat <- data()
      matfilt <- mat[,1:154]
      
      tcounts <- t(matfilt) %>%
        base::merge(bulkseq_colData, ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
      
      ## add gene symbols for goi terms
      tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol
      
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol) %>% dplyr::summarise(expression=median(expression))
      tcounts_med <- tcounts_med[!tcounts_med$Condition %in% "Ipsi", ] #remove all injured samples
      
      g <- ggplot(tcounts_med, aes(x=interaction(Population), y = symbol)) #group by population
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
      
      ggplotly(g) 
    })
    
    output$bulkseq_lines <- renderPlotly({
        
      req(data());
      
      mat <- data()
      matfilt <- mat[,1:154]
      
      tcounts <- t(matfilt) %>%
        base::merge(bulkseq_colData, ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
        
        ## add gene symbols for goi terms
        tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol
        
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, Timepoint, symbol) %>% dplyr::summarise(expression=median(expression))
        
        g <- ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint)) ) 
        g <- g + scale_colour_viridis(discrete=TRUE, end = .80)
        g <- g + geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol))
        #g <- g + facet_wrap(~symbol, ncol=3)
        g <- g + theme_bw() + theme(
            panel.grid = element_blank(),
            #axis.title = element_blank(),
            axis.text.x = element_text(size=10, angle = 45, hjust= 1),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank())

        ggplotly(g)
    })
    
    
    output$bulkseq_lines_subtype <- renderPlotly({
      
      req(data());
      
      mat <- data()
      matfilt <- mat[,1:154]
      
      tcounts <- t(matfilt) %>%
        base::merge(bulkseq_colData, ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
      
      ## add gene symbols for goi terms
      tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol
      
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Timepoint, Population, symbol) %>% dplyr::summarise(expression=median(expression))
      
      g <- ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint)) ) 
      g <- g + scale_colour_viridis(discrete=TRUE, end = .80)
      g <- g + geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol))
      g <- g + facet_wrap(~Population, ncol=5)
      g <- g + theme_bw() + theme(
        panel.grid = element_blank(),
        #axis.title = element_blank(),
        axis.text.x = element_text(size=10, angle = 45, hjust= 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
      
      ggplotly(g)
    })
    
    
    
    
    ### leaflet map for contact details
    output$myMap <- renderLeaflet({
        m <- leaflet() %>% addTiles()
        m <- m %>% setView( -1.238233, 51.756192, zoom = 13)
        m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
    })
    
    #volcanoServer("volcano_module")
    
    
})
