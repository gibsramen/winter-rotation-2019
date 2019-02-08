library(CancerSubtypes)

# ---- LOADING ----

setwd('~/projects/rotation_2019/cancer_subtypes')
print('Loading data frames...')
cib_df <- read.csv('../data/cib_data_sid_consolidated.csv')
otu_df <- read.csv('../data/snmData_Cib_filt.csv')
print('Data frames loaded!')

# ---- FILTERING ----

BRCA_cib_df <- cib_df[cib_df$CancerType == 'BRCA',]
BRCA_cib_df <- BRCA_cib_df[order(BRCA_cib_df$feature.id),]

BRCA_otu_df <- otu_df[otu_df$Unnamed..0 %in% BRCA_cib_df$feature.id,]
BRCA_otu_df <- BRCA_otu_df[order(BRCA_otu_df$Unnamed..0),]

#BRCA_cib_df <- BRCA_cib_df[,!(colnames(BRCA_cib_df) %in% c('feature.id', 'CancerType'))]
#BRCA_otu_df <- BRCA_otu_df[,!(colnames(BRCA_otu_df) %in% c('Unnamed..0'))]

rownames(BRCA_cib_df) <- NULL
rownames(BRCA_otu_df) <- NULL

# ---- TRANSPOSE ----

t.BRCA_cib_df <- t(BRCA_cib_df)
colnames(t.BRCA_cib_df) <- t.BRCA_cib_df[1,]
t.BRCA_cib_df <- t.BRCA_cib_df[-c(1,2),]
class(t.BRCA_cib_df) <- 'numeric'

t.BRCA_otu_df <- t(BRCA_otu_df)
colnames(t.BRCA_otu_df) <- t.BRCA_otu_df[1,]
t.BRCA_otu_df <- t.BRCA_otu_df[-c(1),]
class(t.BRCA_otu_df) <- 'numeric'

# ---- SNF-CC ----

BRCA <- list(OTU=t.BRCA_otu_df, Cib=t.BRCA_cib_df)
result <- ExecuteSNF.CC(datasets = BRCA, clusterNum = 3)
print('Finished SNF & CC!')

# ---- PLOT ----

sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
plot(sil, border=NA, col=2:4)
