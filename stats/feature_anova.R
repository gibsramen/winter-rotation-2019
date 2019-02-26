library(dplyr)
library(reshape2)

# take in group file and cibersort csv as inputs
# return ANOVA of all microbes/cibersort by group
# probably only return significant values

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide arguments!')
}

group_file <- args[1]
cib_file <- args[2]
out_file <- args[3]

print('Loading data frames...')
#otu_df <- read.csv('../data/snmData_Cib_filt.csv')
cib_df <- read.csv(cib_file)
meta_df <- read.csv('../data/metadata_sid_filt.csv')
group_df <- read.csv(group_file, header=F)

t.group_df <- t(group_df)
colnames(t.group_df) <- c('id', 'index', 'group')
t.group_df <- t.group_df[,c('id', 'group')]
t.group_df <- data.frame(t.group_df)
t.group_df <- t.group_df[order(t.group_df$id),]

print('Filtering data frames...')
cib_df <- cib_df[cib_df$feature.id %in% t.group_df$id,]
cib_group_df <- data.frame(cib_df)
cib_group_df <- cib_group_df[order(cib_group_df$feature.id),]
cib_group_df$group <- t.group_df$group
cib_group_df$feature.id <- NULL
cib_group_df$CancerType <- NULL
melt.cib_group_df <- melt(cib_group_df, id.vars = c('group'))

for (cell_type in colnames(cib_group_df)){
	if (cell_type == 'group'){ next }
	print(cell_type)
	cell.cib_group_df <- cib_group_df[,c('group', cell_type)]
	colnames(cell.cib_group_df) <- c('group', 'value')
	aov <- aov(value ~ group, data = cell.cib_group_df)
	print(summary(aov))
	test = summary(aov)
	print(test[[1]]$"Pr(>F)")
	tk <- TukeyHSD(aov)
	print(tk)
}

