##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

#rsconnect::deployApp()
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.15/bioc"))

library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFeedback)
library(plotly)
library(tidyverse)
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("sva")
#library(Tmisc)
#library(calibrate)
library("biomaRt")
#library("IHW")
#library(lattice)
library(plotrix)
library(ggrepel)
library(ggbiplot)
#library(DT)
library(leaflet)
#library(GenomicFeatures)
library(stringr)
library("ComplexHeatmap")
library(ggplot2)
library(viridis)

### Modules

### IU 

shinyUI(fluidPage(
    
  
    #CSS style sheet
    includeCSS("www/style.css"),
    
    shinyFeedback::useShinyFeedback(),
    
    # UI built with shiny dashboard, requires boxes/columns in tabItems 
    dashboardPage(
        dashboardHeader(title="DRG Directory", titleWidth = 225),
        dashboardSidebar(width = 225,
                         sidebarMenu(
                             id = "tabs",
                             
                             HTML(paste0( # oxford logo + ndcn link    
                                 "<br><br>",
                                 "<a href='https://www.ndcn.ox.ac.uk/research/neural-injury-group' target='_blank'> <img style = 'display: block; margin-left: auto; margin-right: auto;' src='oxfordlogo2.png' width = '125'></a>",
                                 "<br><br>")
                             ),
                             
                             # sidebar menu for tabs (pages)    
                             menuItem("Home", tabName = "tabhome", icon = icon("home")),
                             menuItem("Data Tables", tabName = "tabdata", icon = icon("stats", lib = "glyphicon")),
                             menuItem("Dataset details", tabName = "tabsummary", icon = icon("book")),
                             #menuItem("Pathway analysis", tabName = "tabpaths", icon = icon("code-branch")),
                             menuItem("User Guide + Data", tabName = "tabcode", icon = icon("tasks")),
                             menuItem("Contact", tabName = "tabguide", icon = icon("info-circle")),
                             br(),
                             br()
                         ),
                         
                         #Footer (icons + social media links)
                         HTML(paste0(
                             "<table style='clear: both;padding: 0;text-align: center;vertical-align: middle;line-height: normal;
                margin: 30px;position: fixed;bottom: 0px;width: 165px'>", # start table
                
                # icons placed in 6 column table
                "<tr >",
                "<td style='padding: 10px;'></td>",
                "<td style='padding: 5px;'><a href='https://twitter.com/aliibarry' target='_blank'><i class='fab fa-twitter fa-lg'></i></a></td>",
                "<td style='padding: 5px;'><a href='https://github.com/aliibarry' target='_blank'><i class='fab fa-github fa-lg'></i></a></td>",
                "<td style='padding: 5px;'><a href='https://orcid.org/0000-0002-6787-6889' target='_blank'><i class='fab fa-orcid fa-lg'></i></a></td>",
                "<td style='padding: 5px;'><a href='https://www.ndcn.ox.ac.uk/team/allison-barry' target='_blank'><i class='fas fa-brain fa-lg'></i></a></td>",
                "<td style='padding: 10px;'></td>",
                "</tr>",
                
                # second table row, merged columns    
                "<tr>",
                "<script>","var today = new Date();","var yyyy = today.getFullYear();","</script>",
                "<td  colspan='6', style = 'text-align: center;'><small>&copy; - <a href='https://github.com/aliibarry?tab=projects' target='_blank'>github.com/aliibarry</a> - <script>document.write(yyyy);</script></small></td>",
                "</tr>",
                
                "</table>", #end table
                "<br>")
                         )
        ),
        
        # Main panels for output, each TabItem is a different page, same sidebar menu
        dashboardBody(
            tabItems(
                tabItem(tabName="tabhome", 
                        #h4("Home"),
                        fluidRow(
                            
                            # text summary for page
                            box(width=12,
                                status = "primary", 
                                solidHeader = TRUE,
                                #height = 275,
                                h4("Overview"),
                                p("This database provides an interface to explore murine -omics datasets. 
                                  Deep RNA-seq of mouse DRG subpopulations after spare nerve injury (SNI) was performed 
                                  to interrogate subtype-specific and shared injury signatures. 
                                  Specific transgenic details and methodologies will be available in the published reports, 
                                  and are currently available by request. 
                                  All work presented here is currently unpublished, and is the work of Ali Barry, 
                                  Giorgos Baskozos, and Dave Bennett at the University of Oxford (NDCN).
                                Full descriptions and relevant citations can be found in the 'Dataset' tab.  
                                ")
                              )
                            ), #fluidrow
                        
                        fluidRow(
                            column(6, offset = 0,
                                       br(),
                                       h4("Search"),
                                       selectizeInput(
                                         inputId = "geneid", 
                                         label = "", 
                                         multiple = TRUE,
                                         choices = NULL
                                         )
                                       # ,
                                       # br(),
                                       # 
                                       # h4("Search by file (mouse ensembl id):"), 
                                       # fileInput("user_file", 
                                       #           label = NULL,
                                       #           accept = c(
                                       #               'text/csv',
                                       #               'text/comma-separated-values',
                                       #               'text/tab-separated-values',
                                       #               'text/plain',
                                       #               '.csv',
                                       #               '.tsv')
                                       # )
                                ),
                            column(6,
                                   br(),
                                   img(src = "schematic.png", height = 150, width = 400, align = "right")
                              
                              
                            )
                            ),

                        
                        fluidRow(
                            column(12, 
                                   hr(),
                                   h4("Results"),
                                   actionLink("link_to_tables", "Data Tables"),
                                   
                                   p("All results are plotted as median vst-transformed count data. 
                                     Search result data is available for download in in the 'Data Table' tab. 
                                     Interactive queries for hypothesis testing are forthcoming, 
                                     but the csv files are available on github.")
                                   ),
                            
                            column(width = 6,

                                h4("Naive"),
                                plotlyOutput("bulkseq_dots")
                                ),
                            column(width = 6,
                                h4("SNI"),
                                plotlyOutput("bulkseq_lines")
                                )

                            ), #fluidRow
                        fluidRow(
                          column(width = 12,
                                 hr(),
                                 h4("Subtype Plots"),
                                 plotlyOutput("bulkseq_lines_subtype")
                                )
                        )

                ), # tabItem HOME
                
                
                tabItem(tabName="tabdata",
                        fluidRow(
                          column(12, h4("Results")),
                          
                          column(width = 12, 
                                 includeMarkdown("datatable_notes.Rmd")),
                          br(),
                          
                          column(width = 12, 
                                #h4("Table"),
                                downloadButton("downloadData", "Download"),
                                br(),
                                br(),
                                DT::dataTableOutput("goi_table")
                                #DT::dataTableOutput("geneidtable")
                                )
                        )
                          
                        #   column(width=12, 
                        #          downloadButton("downloadData", "Download"))
                        # )
                        
                ), # tabItem one
                
                
                tabItem(tabName="tabsummary",
                        h4("Available datasets"),
                        includeMarkdown("datasetsummary.md")
                ), 
                
                # 
                # ## Pathway analyses for various datasets. GO-term searches, etc.
                # tabItem(tabName="tabpaths",
                #         h4("Pathway and network analyses"),
                #         p("beta")
                #         
                # ), # tabItem paths
                
                ## Supply simple links for each paper + supplementary repository
                tabItem(tabName="tabcode",
                        h4("Data Access"),
                        includeMarkdown('codedata.md'),
                        br(),
                        
                        # HTML links to funders, embedded in a table for alignment
                        HTML(paste0(
                            "<table style='margin-left:auto; margin-right:auto; bottom=0;'>",
                            
                            "<tr>",
                            "<td><a href='https://www.gtc.ox.ac.uk/' target='_blank'> <img style= 'display: center;' src='gtclogo.png', width = '100'></a></td>",
                            "<td><a href='https://wellcome.org/' target='_blank'> <img style = 'display: center;' src='wt-logo.svg'. width = '100'></a></td>",
                            "</tr>",
                            
                            "</table>"
                        )
                        )
                        
                ), #tabitem code
                
                tabItem(tabName="tabguide",
                        h4("Contact"),
                        includeMarkdown('userguide.md'),
                        
                        # Lab location map, in a column solely for aesthetic
                        column(10, offset = 1,
                               leafletOutput('myMap', width = "100%", height = 350)
                        )
                
                ) # tabItem help guide
            ) # tabItems (all) 
        )  #mainpanel dashboardbody close
    ) # dashboad page close
) #fluid page
)#ui
