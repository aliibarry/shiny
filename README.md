# DRG directory

Shiny (R) app to interactively host -omics data.

Hosted at:
https://livedataoxford.shinyapps.io/drg-directory/

***
The userguide is associated with the paper Barry et al  


Allison M Barry, Na Zhao, Xun Yang, David L Bennett, Georgios Baskozos  
‘Deep RNA-seq of male and female murine sensory neuron subtypes after nerve injury’. 

Data currently availale from Barry 2022 (PhD Thesis), [GSE216444](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216444) with PMID to follow.
  
***
  
USER GUIDE: 
1. The data is in .Rdata format. To load user's own data, substitute the directory in "data_dir" variable at the top of server.R file. 
2. The .Rdata should include three files: 
    1) a count matrix, here bulkseq_mat 
    2) a processed differential gene data file that contain the following columns: ('symbol', 'loadfold2change', 'padj') with gene ID as row names. 
    3) a colData, (here bulkseq_colData), describing the Timepoint, Population, Sex, Condition of each sample (or associated factors).  
3. If user's datafile have different names, use the find-replace method to change the name of variables in the file into user's own filenames 
4. Variable names: 
      ** tcounts_med is a data matrix containing median expression of each gene
      ** population_labels, subpopulation_labels and sex labels are full names of abbreviations used; can be substituted according to user's data 
      ** subpopulations is a list of csv file names storing processed datasets for each subpopulation. Can be adjusted based on user's data  
5. There are server modules for each graph type 
      1) df(data, type): a function for processing count data; converts count matrix into a transpose version i.e. row names = sample, and colnames = genes  
      2) plotline_server(id, df, sex): a server module for plotting line graphs; used to show interaction between timepoints and injury. "df" is the count data for selected genes.
              ** if user data does not contain 'Timepoint' variable, add a new column named Timepoint in user's colData with values 'na' for each sample. This will 
              prevent errors when rouping samples by Sex Timepoint; Also, delete "linetype=Timepoint" in the ggplot object. 
      3) plotdot_server(id, df, sex): a server module for plotting Gene expression dot plots showing
              ** plotting gene expression for selected genes.  
              ** if user's data does not contain 'Sex' variable, add a new column named "Sex" in user's colData with values 'mixed' for each sample. This will 
              prevent errors when grouping samples by Sex. Delete "facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels))" in the ggplot
              ** "df" is the count data for selected genes.
              ** Dots are coloured by gene expression.
              ** Dot sizes reflect differences in more highly expressed genes.
      4) plotsubtype_server(id, df, sex): Gene expression dot plots showing differences across population 
              ** if user's data does not contain 'Sex' variable, add a new column named "Sex" in user's colData with values 'mixed' for each sample. This will 
              prevent errors when grouping samples by Sex. Delete "facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels))" in the ggplot
              ** if user's data does not contain 'Population' variable, this plot is useful to be included 
      5) contrast_table_server(id, df): A module for displaying tables showing the significance of load2foldchange. "df" is the processed differential gene data file from DESeq2. 
      6) goi_table_server(id, df): A module for displaying the individual count for selected genes of all samples. "df" is the count data for selected genes. 
      7) deg_plot_server(id, df): A module for plotting the significance of load2foldchange for each gene. "df" is the count data for selected genes.
      8) volcano_plot(df, ils): A function for plotting volcano plots. This function is called in the volcano plot server. "ils" is the list of selected genes.  
              ** Volcano plots are generated per subtype, showing large fold changes in injury(or other conditions)-related genes. 
      9) volcano_plot_server(id, df, ils): A module for plotting volcano plots 
              ** processed differential gene data files are loaded. A volcano plot was plotted for each subtype to show fold changes in injury(or other conditions)-related genes. 
              ** volcano plots are generated using the external function "volcano_plot()"

Additional Information:

   1) If users would like to delete existing plots, they can directly delete server modules. 
   2) If users would like to add more plots, they can create new server modules. 
