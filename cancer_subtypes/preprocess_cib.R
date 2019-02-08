cib_df <- read.csv('../data/cib_data_sid.tsv', sep='\t')

cib_df$B.cells <- cib_df$B.cells.memory + cib_df$B.cells.naive
cib_df$T.cells.CD4 <- cib_df$T.cells.CD4.memory.activated + cib_df$T.cells.CD4.memory.resting +
  cib_df$T.cells.CD4.naive + cib_df$T.cells.follicular.helper +
  cib_df$T.cells.regulatory..Tregs.
cib_df$NK.cells <- cib_df$NK.cells.activated + cib_df$NK.cells.resting
cib_df$Macrophages <- cib_df$Macrophages.M0 + cib_df$Macrophages.M1 + cib_df$Macrophages.M2 + 
  cib_df$Monocytes
cib_df$Mast.cells <- cib_df$Mast.cells.activated + cib_df$Mast.cells.resting
cib_df$Dendritic.cells <- cib_df$Dendritic.cells.activated + cib_df$Dendritic.cells.resting

drops <- c('B.cells.memory', 'B.cells.naive', 'T.cells.CD4.memory.activated',
           'T.cells.CD4.memory.resting', 'T.cells.CD4.naive', 'NK.cells.activated',
           'NK.cells.resting', 'Macrophages.M0', 'Macrophages.M1', 'Macrophages.M2',
           'Monocytes', 'Mast.cells.activated', 'Mast.cells.resting',
           'Dendritic.cells.activated', 'Dendritic.cells.resting',
           'T.cells.regulatory..Tregs.', 'T.cells.follicular.helper')

cib_df_2 <- cib_df[,!colnames(cib_df) %in% drops]
write.csv(cib_df_2, '~/data/TCGA_immune_microbiome/cib_data_sid_consolidated.csv', row.names=FALSE)
print('Written!')
