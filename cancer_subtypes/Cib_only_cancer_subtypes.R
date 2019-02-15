library(CancerSubtypes)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# input cancer type by which to cluster
cancer_type <- args[1]
print(cancer_type)

# ---- LOADING ----

setwd('~/projects/rotation_2019/cancer_subtypes')
print('Loading data frames...')
cib_df <- read.csv('../data/cib_data_sid_consolidated_clr.csv')
print('Data frames loaded!')

# ---- FILTERING ----

cancer_cib_df <- cib_df[cib_df$CancerType == cancer_type,]
if (dim(cancer_cib_df)[1] < 100){
	stop('Cancer has fewer than 100 samples!')
}
cancer_cib_df <- cancer_cib_df[order(cancer_cib_df$feature.id),]

rownames(cancer_cib_df) <- NULL
print('Finished filtering data frames!')

# ---- TRANSPOSE ----
# SNF requires rows as features and columns as samples (opposite csv)
# first row in transposed df is sample
# in cib_df, second row is cancer type
# by default transpose seems to coerce to characters

t.cancer_cib_df <- t(cancer_cib_df)
colnames(t.cancer_cib_df) <- t.cancer_cib_df[1,]
t.cancer_cib_df <- t.cancer_cib_df[-c(1,2),] 
class(t.cancer_cib_df) <- 'numeric'
t.norm.cancer_cib_df <- data.normalization(t.cancer_cib_df,
	type="feature_zscore", log2=FALSE)

# ---- SNF-CC ----

avg_sil_widths = c()
max_num <- 10
for (i in c(2:max_num)){
	print(sprintf('Testing %i clusters...', i))
	out_dir = paste('cib_only_cluster_results', cancer_type, 
					sprintf('%i_clusters', i), sep='/')
	dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
	# using list of 1 dataset - normalizes to z-score?
	result <- ExecuteCC(clusterNum = i, d=t.norm.cancer_cib_df,
		title=out_dir, plot='png')
	print('Finished SNF & CC!')
	group <- result$group
	group.df <- data.frame(rbind(cancer_cib_df$feature.id, group))
	write.csv(group.df, paste(out_dir, 'group.csv', sep='/'), row.names=FALSE)

	sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
	png(paste(out_dir, 'silhouette.png', sep='/'))
	plot(sil, border=NA, col=2:(i+1))
	dev.off()

	avg_sil_widths[i-1] = mean(sil[,3])
}

sil_out_file = paste('cib_only_cluster_results',
	cancer_type, 'silhouette_widths.tsv',  sep='/')
sil.df <- data.frame(cbind(c(2:max_num), avg_sil_widths))
colnames(sil.df) <- c('Num_Clusters', 'Avg_Sil_Width')
write.table(sil.df, sil_out_file, row.names = FALSE, sep='\t')
