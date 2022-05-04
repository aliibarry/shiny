
save(file = paste0(PATH_results, "/drg-directory.RData"), 
     bulkseq_colData, bulkseq_mat, gene_data, PATH_results
)

PATH_results = "~/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/amb/Database/drg-directory"

bulkseq_colData <- colData(bulkseq_dds)
