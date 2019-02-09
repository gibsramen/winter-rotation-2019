library(CancerSubtypes)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

cancer_type <- args[1]

# ---- LOADING ----

setwd('~/projects/rotation_2019/cancer_subtypes')
print('Loading data frames...')
cib_df <- read.csv('../data/cib_data_sid_consolidated_clr.csv')
otu_df <- read.csv('../data/snmData_Cib_filt.csv')
meta_df <- read.csv('../data/metadata_sid_filt.csv')
print('Data frames loaded!')

# ---- FILTERING ----

cancer_cib_df <- cib_df[cib_df$CancerType == cancer_type,]
if (dim(cancer_cib_df)[1] < 100){
	stop('Cancer has fewer than 100 samples!')
}
cancer_cib_df <- cancer_cib_df[order(cancer_cib_df$feature.id),]

cancer_otu_df <- otu_df[otu_df$Unnamed..0 %in% cancer_cib_df$feature.id,]
cancer_otu_df <- cancer_otu_df[order(cancer_otu_df$Unnamed..0),]

cancer_meta_df <- meta_df[meta_df$X %in% cancer_cib_df$feature.id,]
cancer_meta_df <- cancer_meta_df[order(cancer_meta_df$X),]

cancer_meta_df <- cancer_meta_df[,colnames(cancer_meta_df) %in% c('X', 'race', 'age_at_diagnosis', 'gender')]

cancer_race_df <- cancer_meta_df[,colnames(cancer_meta_df) %in% c('X', 'race')]
cancer_race_df <- cancer_race_df[order(cancer_race_df$X),]

cancer_age_df <- cancer_meta_df[,colnames(cancer_meta_df) %in% c('X', 'age_at_diagnosis')]
cancer_age_df <- cancer_age_df[order(cancer_age_df$X),]

cancer_gender_df <- cancer_meta_df[,colnames(cancer_meta_df) %in% c('X', 'gender')]
cancer_gender_df <- cancer_gender_df[order(cancer_gender_df$X),]

rownames(cancer_cib_df) <- NULL
rownames(cancer_otu_df) <- NULL
rownames(cancer_meta_df) <- NULL
rownames(cancer_race_df) <- NULL
rownames(cancer_gender_df) <- NULL
rownames(cancer_age_df) <- NULL
print('Finished filtering data frames!')

# ---- TRANSPOSE ----

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
result <- ExecuteSNF.CC(datasets = cancer_data, clusterNum = 3, 
						title = cancer_type, plot = 'png')
group <- result$group
group.df <- data.frame(rbind(cancer_cib_df$feature.id, group))
write.csv(group.df, paste(cancer_type, 'group.csv', sep='/'), row.names=FALSE)
print('Finished SNF & CC!')

# ---- PLOT ----

sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
png(paste(cancer_type, 'silhouette.png', sep='/'))
plot(sil, border=NA, col=2:4)
dev.off()
