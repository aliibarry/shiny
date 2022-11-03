##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################
library(profvis)
# global variables 
data_dir = "/Users/lynnezhao/Desktop/drg/drg-directory.RData"
load(data_dir)


# mat = bulkseq_mat # count data 
population_labels = c("TBAC" = "A\u03b4-LTMR + A\u03b2-RA-LTMR",
                      "CRTH" = "C-LTMR",
                      "MRTD" = "NP",
                      "CGRT" = "PEP",
                      "TDNV" = "Nociceptors")
sexlabels = c("F" = "Female", "M" = "Male")
subpopulations = c(
  "Nociceptors 3D" = "TDNV_3D.csv",
  "Nociceptors 4W" = "TDNV_4W.csv",
  "PEP 3D" = "CGRT_3D.csv",
  "PEP 4W" = "CGRT_4W.csv", 
  "NP 3D" = "MRTD_3D.csv",
  "NP 4W" = "MRTD_4W.csv",
  "C-LTMR 3D" = "CRTH_3D.csv",
  "C-LTMR 4W" = "CRTH_4W.csv",
  "Ad- AB-RA LTMRs 3D" = "TBAC_3D.csv",
  "Ad- AB-RA LTMRs 4W" = "TBAC_4W.csv"
)

subpopulation_labels = c(
  "TDNV_3D" = "Nociceptors 3D",
  "TDNV_4W" = "Nociceptors 4W",
  "CGRT_3D" = "PEP 3D",
  "CGRT_4W" = "PEP 4W", 
  "MRTD_3D" = "NP 3D",
  "MRTD_4W" = "NP 4W",
  "CRTH_3D" = "C-LTMR 3D",
  "CRTH_4W" = "C-LTMR 4W",
  "TBAC_3D" = "Ad- AB-RA LTMRs 3D",
  "TBAC_4W" = "Ad- AB-RA LTMRs 4W"
)

# a function that parse data into different dataframes based on the graph type 
df <- function(data, type) {
  matfilt <- data()[,1:154]
  tcounts <- t(matfilt) %>%
    base::merge(bulkseq_colData, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
  tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
  return(tcounts)
}

# a server module for ploting line plots 
plotline_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    if (sex == 'Both') {
      tcounts_med <- df %>% dplyr::group_by(Condition, Timepoint, symbol) %>% 
        dplyr::summarise(expression=median(expression))
      output$bulkseq_lines <- renderPlotly({
        ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) + 
          scale_colour_viridis(discrete=TRUE, end = .80) + 
          geom_line(aes(color=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_bw() + 
          theme(panel.grid = element_blank(),
            axis.title = element_text(family = "Gadugi", size=12),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(), legend.justification = c(0,0.3)) + ylab("Expression") +
          guides(col=FALSE, linetype=guide_legend(label.hjust=0.5, ncol=8)) 
      })
    }
    
    if (sex =='Separate'){
      tcounts_med <- df %>% dplyr::group_by(Condition, Timepoint, symbol, Sex) %>% 
        dplyr::summarise(expression=median(expression))
      output$bulkseq_lines <- renderPlotly({
        ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint)) ) + 
          scale_colour_viridis(discrete=TRUE, end = .80) + 
          facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels))+
          geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_bw() + 
          theme(
            panel.grid = element_blank(),
            axis.title = element_text(family = "Gadugi", size=12),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) + ylab("Expression")+
          labs(col="")
      })
    }
  })
}


# a server module for ploting dot plots 
plotdot_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    if (sex == 'Both') {
      tcounts_med <- df %>% dplyr::group_by(Condition, Population, symbol) %>% 
        dplyr::summarise(expression=median(expression))
      tcounts_med <- tcounts_med[!tcounts_med$Condition %in% "Ipsi", ] #remove all injured samples
      output$bulkseq_dots <- renderPlotly({
        ggplot(data = tcounts_med, aes(x= Population, y = symbol)) + 
          scale_colour_viridis_c(option = "magma", end = .90) + 
          geom_point(aes(col=expression, size=expression)) + theme_bw() + theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) + scale_x_discrete(labels=population_labels)
      })
    }
    if (sex == "Separate") {
      tcounts_med <- df %>% dplyr::group_by(Condition, Population, symbol, Sex) %>% 
        dplyr::summarise(expression=median(expression))
      tcounts_med <- tcounts_med[!tcounts_med$Condition %in% "Ipsi", ] #remove all injured samples
      output$bulkseq_dots <- renderPlotly({
        ggplot(data = tcounts_med, aes(x= Population, y = symbol)) + 
          scale_colour_viridis_c(option = "magma", end = .90) + 
          facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels))+
          geom_point(aes(col=expression, size=expression)) + theme_bw() + theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) + scale_x_discrete(labels=population_labels)
      })
    }
      })
}

# a server module for ploting subtype plots; 
plotsubtype_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    if (sex == 'Both') {
      tcounts_med <- df %>% dplyr::group_by(Condition, Timepoint, Population, symbol) %>% 
        dplyr::summarise(expression=median(expression))
      output$bulkseq_lines_subtype <- renderPlotly({
        ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +  
          scale_colour_viridis(discrete=TRUE, end = .80) + 
          geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + 
          # facet_grid(symbol~Population, labeller=labeller(Population=population_labels)) +
          facet_wrap(~Population, ncol=5, labeller=labeller(Population=population_labels))+ 
          theme_bw() + theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.title = element_text(family = "Gadugi", size=12),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) + ylab("Expression") + labs(col="")
      })
    }
    if (sex=='Separate'){
      tcounts_med <- df %>% dplyr::group_by(Condition, Timepoint, Population, symbol, Sex) %>% 
        dplyr::summarise(expression=median(expression))
      output$bulkseq_lines_subtype <- renderPlotly({
        ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +  
          scale_colour_viridis(discrete=TRUE, end = .80) + 
          geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + 
          # facet_grid(symbol~Population, labeller=labeller(Population=population_labels)) +
          facet_grid(Sex~Population, labeller=labeller(Population=population_labels, Sex=sexlabels))+ 
          theme_bw() + theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
            axis.text.y = element_text(family = "Gadugi", size=10),
            axis.title = element_text(family = "Gadugi", size=12),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()) + ylab("Expression") + labs(col="")
      })
    }
   
  })
}

# a server for rendering contrast tables 
contrast_table_server <- function(id, df) {
  moduleServer(id,function(input,output,session) {
    output$contrast_table <- DT::renderDataTable({
      req(df())
      DT::datatable(
        df(),
        width = 12,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
  })
}

# a server for rendering goi tables 
goi_table_server <- function(id, df) {
  moduleServer(id, function(input, output, session){
    output$goi_table <- DT::renderDataTable({
      req(df())
      datatable <- df()
      datatable$geneID <- rownames(datatable)
      DT::datatable(
        datatable,rownames=datatable$symbol,
        style="default",
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })
  })
}

# plot deg servers
deg_plot_server <- function(id, df) {
  moduleServer(id, function(input, output, session){
    output$deg_plot <- renderPlotly({
      ggplot(data = df(), aes(x=interaction(Population), y = symbol, 
                              text = paste('FDR: ',FDR))) + 
        # facet_grid(species~dataset, scale="free", space="free") +
        scale_colour_viridis_c(option = "viridis", end = .90) +
        geom_point(aes(col=LFC, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
        theme_bw() + theme(
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(family = "Gadugi", size=8, angle = 45, hjust= 1),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) + 
        labs(shape = "", size = "") + scale_x_discrete(labels=subpopulation_labels)
    })

  })
}

# functions for ploting volcano plots + highlighting selected genes
volcano_plot <- function(df, ils) {
  
  il_genes <- df %>% filter(symbol %in% ils) 
  
  ggplot(data = df, aes(x = LFC, y = log10fdr)) + 
    geom_point(colour = "grey", alpha = 0.5) +
    geom_point(data = il_genes, # New layer containing data subset il_genes       
               size = 2,
               shape = 21,
               fill = "firebrick",
               colour = "black") + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    geom_label_repel(data = il_genes,   
                     aes(label = symbol),
                     force = 2,
                     nudge_y = 1) + 
    scale_colour_manual(values = cols) + 
    scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                       limits = c(-10, 10)) + theme_bw() + theme(
                         panel.grid = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank(), 
                         axis.text.x = element_text(size=12)) +
    labs(y="-log10(FDR)")
}

volcano_plot_server <- function(id, ils) {
  moduleServer(id, function(input, output, session){
    
    output$volcano1 <- renderPlot({
      res = fread("TDNV_3D.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano2 <- renderPlot({
      res = fread("TDNV_4W.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano3 <- renderPlot({
      res = fread("CGRT_3D.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano4 <- renderPlot({
      res = fread("CGRT_4W.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano5 <- renderPlot({
      res = fread("MRTD_3D.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano6 <- renderPlot({
      res = fread("MRTD_4W.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano7 <- renderPlot({
      res = fread("CRTH_3D.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano8 <- renderPlot({
      res = fread("CGRT_4W.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano9 <- renderPlot({
      res = fread("TBAC_3D.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
    output$volcano10 <- renderPlot({
      res = fread("TBAC_4W.csv")
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      mutated_ress = mutate(res, log10fdr=-log10(FDR))
      volcano_plot(mutated_ress, ils)
    })
    
  })
}

shinyServer(function(input, output, session) {
  load(data_dir)
  updateSelectizeInput(session, 
                       inputId = "geneid", 
                       label = "Search Genes:",
                       choices = bulkseq_mat[,155], 
                       server = TRUE,
                       selected = "Atf3"
                       #options = list(placeholder = 'select a gene', 'plugins' = list('remove_button'))
                       ) # select genes 
  
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
    
    observeEvent(input$link_to_home, {
      newvalue <- "tabhome"
      updateTabItems(session, "tabs", newvalue)
    })
    
    observeEvent(input$link_to_home2, {
      newvalue <- "tabhome"
      updateTabItems(session, "tabs", newvalue)
    })
    
    # data preprocessing for deg plots 
    deg_df <- data.frame()
    for (pop in subpopulations) {
      res <- fread(pop)
      colnames(res) <- c('symbol', 'LFC', 'FDR', 'ID')
      res$Population <- rep(substring(pop, 1, 7), nrow(res)) # add the population label
      deg_df = bind_rows(deg_df, res)
    }
    mutateddf = mutate(deg_df, sig=ifelse(deg_df$FDR<0.05, "SIG", "NS"))
   
    # ploting volcano plots
    observeEvent(input$volc, {
      # plot volcano plots only when clicking the button
      volcano_plot_server("volcano", input$geneid)
    })
    
    # for updating selected genes and plotting graphs, only after clicking 'plot graphs'
    observeEvent(input$load, {
      data <- reactive({
        bulkseq_mat[bulkseq_mat[,155] %in% input$geneid,]
      }) %>% bindCache(input$geneid)
      plotdot_server("dot",df(reactive({data()}), "dot"), input$sex) #dotplot
      plotline_server("line", df(reactive({data()}), "line"), input$sex) # lineplot
      plotsubtype_server("lines_subtype", df(reactive({data()}), "lines_subtype"), input$sex) 
      
      # for differential gene analysis plot
      deg <- reactive({
        df <- mutateddf %>% filter(symbol %in% input$geneid)
        return(df)
      }) %>% bindCache(input$geneid)
      
      deg_plot_server("deg_plot", reactive({deg()}))
      goi_table_server("goi_table", reactive({data()}))
      
    })
    
    # display table in 'hypothesis testing' page 
    usercontrast <- reactive({
      req(input$contrast)
      res <- fread(input$contrast)
      return(res)
    }) %>% bindCache(input$contrast)
  
    contrast_table_server("contrast_table",reactive({usercontrast()})) # a contrast table server 

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

    ### leaflet map for contact details
    output$myMap <- renderLeaflet({
        m <- leaflet() %>% addTiles()
        m <- m %>% setView( -1.238233, 51.756192, zoom = 13)
        m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
    })
    
    #volcanoServer("volcano_module")
})

# profiling 
# profvis({
  # shiny::runApp('Desktop/drg')
# })
