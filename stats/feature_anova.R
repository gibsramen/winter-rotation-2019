library(dplyr)
library(reshape2)

# take in group file and cibersort csv as inputs
# return ANOVA of all microbes/cibersort by group
# probably only return significant values

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3){
	stop('Please provide arguments!')
}

group_file <- args[1]
cib_file <- args[2]
out_file <- args[3]

print('Loading data frames...')
#otu_df <- read.csv('../data/snmData_Cib_filt.csv')
cib_df <- read.csv(cib_file)
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

p.vals <- c()
for (cell_type in colnames(cib_group_df)){
	if (cell_type == 'group'){ next }
	print(cell_type)
	cell.cib_group_df <- cib_group_df[,c('group', cell_type)]
	colnames(cell.cib_group_df) <- c('group', 'value')
	aov <- aov(value ~ group, data = cell.cib_group_df)
	summ <- summary(aov)
	p.val <- summ[[1]]$"Pr(>F)"[1]
	p.vals <- c(p.vals, p.val)
	tk <- TukeyHSD(aov)
}

adj.p.vals <- p.adjust(p.vals, method="fdr")

#out.df <- rbind(colnames(cib_group_df), adj.p.vals)
cell_types <- colnames(cib_group_df)[!colnames(cib_group_df) %in% c('group')]
out.df <- data.frame(cbind(cell_types, p.vals, adj.p.vals))
colnames(out.df) <- c('Cell.Type', 'p-vals', 'BH-adj-p-vals')

write.csv(out.df, file = out_file, row.names=F)
