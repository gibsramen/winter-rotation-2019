library(CancerSubtypes)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# input cancer type by which to cluster
cancer_type <- args[1]
cib_file <- args[2]
out_folder <- args[3]
print(cancer_type)

# ---- LOADING ----

setwd('~/projects/rotation_2019/cancer_subtypes')
print('Loading data frames...')

cib_df <- read.csv(cib_file)
otu_df <- read.csv('../data/snmData_Cib_filt.csv')
meta_df <- read.csv('../data/metadata_sid_filt.csv')
days_df <- read.table('../clinical/follow_up_counts.txt', sep='\t')
colnames(days_df) <- c('barcode', 'status', 'days')
days_df$status <- as.character(days_df$status)
print('Data frames loaded!')

# ---- FILTERING ----

cancer_cib_df <- cib_df[cib_df$CancerType == cancer_type,]
if (dim(cancer_cib_df)[1] < 100){
	stop('Cancer has fewer than 100 samples!')
}
cancer_cib_df <- cancer_cib_df[order(cancer_cib_df$feature.id),]

cancer_otu_df <- otu_df[otu_df$Unnamed..0 %in% cancer_cib_df$feature.id,]
cancer_otu_df <- cancer_otu_df[order(cancer_otu_df$Unnamed..0),]

surv_df <- meta_df[,c('X', 'investigation', 'case_id')]
surv_df$investigation <- sub("TCGA-", '', surv_df$investigation)
cancer_surv_df <- surv_df[surv_df$X %in% cancer_cib_df$feature.id,]

status_vector <- c()
days_vector <- c()

days_df$status <- as.character(days_df$status)

for (case in cancer_surv_df$case_id){
	row <- days_df[case == days_df$barcode,]
	status_vector <- c(status_vector, row$status)
	days_vector <- c(days_vector, row$days)
}

status_vector <- ifelse(status_vector == 'Dead', 1, 0)
cancer_surv_df$status <- status_vector
cancer_surv_df$time <- days_vector
cancer_surv_df <- cancer_surv_df[order(cancer_surv_df$X),]

rownames(cancer_cib_df) <- NULL
rownames(cancer_otu_df) <- NULL
rownames(cancer_surv_df) <- NULL
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

t.cancer_otu_df <- t(cancer_otu_df)
colnames(t.cancer_otu_df) <- t.cancer_otu_df[1,]
t.cancer_otu_df <- t.cancer_otu_df[-c(1),]
class(t.cancer_otu_df) <- 'numeric'
print('Finished transposing data frames!')

# ---- SNF-CC ----

cancer_data <- list(OTU=t.cancer_otu_df, Cib=t.cancer_cib_df)

avg_sil_widths = c()
p_vals <- c()
max_num <- 10
for (i in c(2:max_num)){
	print(sprintf('Testing %i clusters...', i))
	#out_dir = paste('cluster_results', cancer_type, 
					#sprintf('%i_clusters', i), sep='/')
	out_dir <- paste(out_folder, cancer_type, sprintf('%i_clusters', i),
		sep = '/')
	dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
	result <- ExecuteSNF.CC(datasets = cancer_data, clusterNum = i,
							title = out_dir, plot = 'png')
	print('Finished SNF & CC!')
	group <- result$group
	group.df <- data.frame(rbind(cancer_cib_df$feature.id, group))
	write.csv(group.df, paste(out_dir, 'group.csv', sep='/'), row.names=FALSE)

	# need to save distance matrix
	dist.matrix <- result$distanceMatrix
	write.csv(dist.matrix, paste(out_dir, 'dist.csv', sep='/'), row.names=FALSE)

	p_1 <- survAnalysis(mainTitle = 'test', cancer_surv_df$time,
						cancer_surv_df$status, group,
						distanceMatrix = dist.matrix,
						similarity = TRUE)
	p_vals <- c(p_vals, p_1)

	sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
	png(paste(out_dir, 'silhouette.png', sep='/'))
	plot(sil, border=NA, col=2:(i+1))
	dev.off()

	avg_sil_widths[i-1] = mean(sil[,3])
}

sil_out_file = paste(out_folder, cancer_type, 'silhouette_widths.tsv', sep='/')
sil.df <- data.frame(cbind(c(2:max_num), avg_sil_widths))
sil.df$p.vals <- p_vals
colnames(sil.df) <- c('Num_Clusters', 'Avg_Sil_Width', 'Surv_P')
write.table(sil.df, sil_out_file, row.names = FALSE, sep='\t')
