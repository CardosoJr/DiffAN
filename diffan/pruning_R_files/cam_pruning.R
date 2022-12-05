print("INICIO DO SCRIPT")
library(mgcv)

source("DiffAN/diffan/pruning_R_files/train_gam.R", chdir=T)
source("DiffAN/diffan/pruning_R_files/selGam.R", chdir=T)
source("DiffAN/diffan/pruning_R_files/pruning.R", chdir=T)



print('SCRIPT R')
dataset <- read.csv(file='{PATH_DATA}', header=FALSE, sep=",")
print(nrow(dataset))
dag <- read.csv(file='{PATH_DAG}', header=FALSE, sep=",")
print(nrow(dag))
set.seed(42)
pruned_dag <- pruning(dataset, dag, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = {CUTOFF}, numBasisFcts = 10), output={VERBOSE})
print(nrow(pruned_dag))

write.csv(as.matrix(pruned_dag), row.names = FALSE, file = '{PATH_RESULTS}')
