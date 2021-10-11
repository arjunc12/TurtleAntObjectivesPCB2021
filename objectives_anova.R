library(rcompanion)
library(FSA)
library(scales)
library(stats)
library(ARTool)
library(dunn.test)
library(PMCMR)
library(conover.test)

args = commandArgs(trailingOnly=TRUE)

percentiles_fname = 'objective_percentiles.csv'

if('--paper' %in% args)
{
    percentiles_fname = 'paper_data/plot_data/objective_percentiles_paper2018-19.csv'
}

reformat_day <- function(day)
{
    afternoon = grepl('afternoon', day)
    day = substr(day, 1, 10)
    if(afternoon)
    {
        day = paste(day, 'pm')
    }
    return(day)
}

objectiveConover <- function(df)
{   
    CT = conover.test(x=df$percentile, g=df$objective, method='bonferroni', kw=FALSE, table=FALSE)
    CT$chi2 = NULL
    CT = as.data.frame(CT)
    CT = CT[,c('comparisons', 'T', 'P', 'P.adjusted')]
    print(CT)
}

df = read.csv(percentiles_fname)
df$day = sapply(df$day, reformat_day)
days = df$day
colonies = df$colony

anova_colonies = c('turtle_hill', 'tejon7', 'T189', 'T500')

anova_objectives = c('total nodes', 'average transition index', 'average length')

df = df[df$colony %in% anova_colonies,]
df = df[df$objective %in% anova_objectives,]

anova_null_model = 'random progression subgraph'
if('--independent' %in% args)
{
    anova_null_model = 'random independent subgraph'   
}
df = df[df$null.model == anova_null_model,]

max_ti = 4
if('--max3' %in% args)
{
    max_ti = 3
}
df = df[df$max.transition.index == max_ti,]

df = df[df$algorithm == 'ants',]


colony_pvals = integer(length(df$day))
colony_objective_pvals = integer(length(df$day))

for(i in 1:nrow(df))
{
    df2 = df[-i,]
    aov = scheirerRayHare(percentile ~ objective + colony, data=df2, verbose=F)
    pvals = aov$p.value
    colony_pval = pvals[1]
    colony_objective_pval = pvals[3]
    colony_pvals[i] = colony_pval
    colony_objective_pvals[i] = colony_objective_pval
}

aov = scheirerRayHare(percentile ~ objective + colony, data=df)
print(aov)


DT1 = dunnTest(percentile ~ objective , data=df, method='bonferroni')
#print(DT1)

objectiveConover(df)

for(colony in anova_colonies)
{
    print(colony)
    colony_df = df[df$colony == colony,]
    
    DT2 = dunnTest(percentile ~ objective, data=colony_df, method='bonferroni')
    #print(DT2)
    
    objectiveConover(colony_df)
}