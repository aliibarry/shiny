##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

#load("C:/Users/allim/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/amb/Database/drg-directory/drg-directory.RData")

#rsconnect::deployApp()
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.15/bioc"))

library(shiny)
library(data.table)
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
#library("ComplexHeatmap")
library(ggplot2)
library(viridis)
library("pacman")
library(profvis)

# global variables 
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
### Modules
plotdot_ui <- function(id, sex) {
  column(width = 6,
         h4("Naive"),
         plotlyOutput(NS(id, "bulkseq_dots"))
  )
}


plotline_ui <- function(id, sex) {
  column(width = 6,
         h4("Injury"),
         plotlyOutput(NS(id, "bulkseq_lines"))
  )
}

plotsubtype_ui <- function(id, sex) {
  column(width = 12,
         h4("Subtype Results"),
         plotlyOutput(NS(id, "bulkseq_lines_subtype"))
  )
}

deg_plot_ui <- function(id) {
  column(width = 6, 
         h4("Differential Gene Analysis"),
         actionLink("link_to_wald", "Differential Gene Analysis Table", icon = icon("stats", lib = "glyphicon")),
         plotlyOutput(NS(id, "deg_plot")))
}

volcano_plot_ui <- function(id) {
  column(width = 6, 
         h4("Ipsilateral vs Contralateral"),
         actionButton("volc", "Plot Volcano Graphs"),
         h5(icon('fas fa-exclamation-triangle'),'WARNING: Plotting may slow server.'),
         tabsetPanel(
           tabPanel("Nociceptors 3D", plotOutput(NS(id, "volcano1"))),
           tabPanel("PEP 3D", plotOutput(NS(id, "volcano3"))),
           tabPanel("NP 3D", plotOutput(NS(id, "volcano5"))),
           tabPanel("C-LTMR 3D", plotOutput(NS(id, "volcano7"))),
           tabPanel("Ad- AB-RA LTMRs 3D", plotOutput(NS(id, "volcano9"))), 
           tabPanel("Nociceptors 4W", plotOutput(NS(id, "volcano2"))),
           tabPanel("PEP 4W", plotOutput(NS(id, "volcano4"))),
           tabPanel("NP 4W", plotOutput(NS(id, "volcano6"))),
           tabPanel("C-LTMR 4W", plotOutput(NS(id, "volcano8"))),
           tabPanel("Ad- AB-RA LTMRs 4W", plotOutput(NS(id, "volcano10")))
         )
  )
}

contrast_table_ui <- function(id) {
  column(12,
         selectInput("contrast", "", 
                     choices = subpopulations,
                     selected = ""),
         DT::dataTableOutput(NS(id, "contrast_table"))
  )
}

goi_table_ui <- function(id) {
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
                   actionLink("link_to_home2", "Back to Plots", icon = icon("home")),
                   br(),
                   br(),
                   DT::dataTableOutput(NS(id,"goi_table"))
            )
          )
  )
}


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
                       menuItem("DEG Analysis", tabName = "tabwald", icon = icon("stats", lib = "glyphicon")),
                       # menuItem("Plots", tabName = "tabplots", icon = icon("chart-bar")),
                       menuItem("Tables", tabName = "tabdata", icon = icon("dna")),
                       menuItem("Dataset summary", tabName = "tabsummary", icon = icon("file-alt")),
                       menuItem("User Guide + Data", tabName = "tabcode", icon = icon("folder-open")),
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
                                  to interrogate neuronal subtype-specific and shared injury signatures. 
                                  Transgenic mice were used to labelling one of five sensory neuron subtypes of interest at 3 days (3D)
                                  and 4 weeks (4W) after injury. 
                                  Populations include: general nociceptors, peptidergic (PEP) and non-peptidergic (NP)
                                  nociceptors, C-LTMRs, and AB-RA + Ad-LTMRs. 
                                  Full descriptions and relevant citations can be found in the 'Dataset' tab.  
                                ")
                  )
                ), #fluidrow
                hr(),
                fluidRow(
                  column(3,offset = 0, 
                         selectizeInput(
                           inputId = "geneid", 
                           label = "Search Genes:", 
                           multiple = TRUE,
                           choices = NULL
                         ), actionButton("load", "Plot Graphs"), br(),br(),
                         actionLink("link_to_tables", "Table: search results", icon = icon("dna"))), 
                  column(3, offset = 0,
                         selectizeInput(
                           inputId = "sex", 
                           label = "Select Sex:", 
                           choices = c('Both', 'Separate'), 
                           selected = 'Both')),
                  br(),br(),
                  br()
                ),
                
                fluidRow(
                  column(width=12, 
                         hr()
                         #h4("Results")
                  ),
                  column(width = 12,
                         p("All results are plotted as median vst-transformed count data. 
                                     Search result data is available for download in in the 'Tables' tab. 
                                     Interactive queries for hypothesis testing with FDR and log2 fold 
                                     change (LFC) details are linked below.")
                  ),
                  plotdot_ui("dot"), 
                  plotline_ui("line"), 
                  plotsubtype_ui("lines_subtype"), 
                  deg_plot_ui("deg_plot"), 
                  volcano_plot_ui("volcano")
                )
                
                
        ), # tabItem HOME
        goi_table_ui("goi_table"), 
        
        tabItem(tabName="tabsummary",
                h4("Available datasets"),
                includeMarkdown("datasetsummary.md"),
                fluidRow(
                  br(),
                  column(12,
                         img(src = "schematic.png", height = 150, width = 400))
                )
        ), 
        
        ## Pathway analyses for various datasets. GO-term searches, etc.
        tabItem(tabName="tabwald",
                h4("DEG Analysis"),
                actionLink("link_to_home", "Back to Plots", icon = icon("home")),
                fluidRow(
                  # deg_plot_ui("deg_plot"),
                  contrast_table_ui("contrast_table")
                )
                
        ),
        
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