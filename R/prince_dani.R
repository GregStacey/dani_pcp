
require(reshape2)
require(ggplot2)
require(readr)
require(PrInCE)
require(dplyr)

first = dplyr::first

# read command line args
this.args = commandArgs(trailingOnly = T)
arg = -1
if (length(this.args) < 1) {
  print("Testing all sets...")
  arg = -1
} else if (length(this.args) == 1 ){
  print("Testing one hyperparameter set...")
  arg = as.integer(as.numeric(this.args))
  ia = (arg %% 7) + 1
  ib = ceiling(arg/7)
}


# load data
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/condition1.csv"
tmp = as.data.frame(read_csv(fn))
tmp = tmp[!grepl("__", tmp$protein.ids), ] # remove REV and CON
#tmp$protein.ids = sapply(sapply(tmp$protein.ids, strsplit, "|", T), "[", 2)
condition1 = (as.matrix(tmp[,3:ncol(tmp)]))
rownames(condition1) = tmp$protein.ids
condition1 = condition1[!is.na(rownames(condition1)),]
condition1 = list(condition1[tmp$replicate==1,], condition1[tmp$replicate==2,], condition1[tmp$replicate==3,])
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/condition2.csv"
tmp = as.data.frame(read_csv(fn))
tmp = tmp[!grepl("__", tmp$protein.ids), ] # remove REV and CON
#tmp$protein.ids = sapply(sapply(tmp$protein.ids, strsplit, "|", T), "[", 2)
condition2 = (as.matrix(tmp[,3:ncol(tmp)]))
rownames(condition2) = tmp$protein.ids
condition2 = condition2[!is.na(rownames(condition2)),]
condition2 = list(condition2[tmp$replicate==1,], condition2[tmp$replicate==2,], condition2[tmp$replicate==3,])
data = list(condition1, condition2)

# load gold_standard
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/allComplexesmapped.txt"
tmp = as.data.frame(read_tsv(fn))
gs.mapped = sapply(tmp$subunits.UniProt.IDs., strsplit, ";")
names(gs.mapped) = tmp$ComplexName
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/allComplexes.txt"
tmp = as.data.frame(read_tsv(fn))
gs = sapply(tmp$`subunits(UniProt IDs)`, strsplit, ";")
names(gs) = tmp$ComplexName
gold_standard = list(gs, gs.mapped)

# replicate combos
reps = list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3), c(1, 2, 3))

# predict
ml.ints = PrInCE(data[[1]][reps[[ia]]], gold_standard[[ib]], verbose = T, classifier="NB")
hl.ints = PrInCE(data[[2]][reps[[ia]]], gold_standard[[ib]], verbose = T, classifier="NB")
sf = paste("../data/interactions_reps",ia,"_gs",ib,".Rda", sep="")
save(ml.ints, hl.ints, file=sf)





################################ Diagnostics

# how many proteins in common b/w gold stanard and data?
unqprots.data = unique(c(rownames(data[[1]][[1]]), rownames(data[[1]][[2]]), rownames(data[[1]][[3]])))
unqprots.gold = unique(unlist(sapply(gold_standard, strsplit, ";")))
n.overlap = length(intersect(unqprots.gold, unqprots.data))
print(paste(length(unqprots.data),"unique proteins in data"))
print(paste(length(unqprots.gold),"unique proteins in gold standard"))
print(paste(n.overlap, " proteins in common"))

# how correlated are true positives? true negatives? random?
uu = 1 # condition
mm = 3 # replicate

# pairwise correlation of all proteins
tmp = data[[uu]][[mm]]
cor.mat = cor(t(tmp), use = "p")
rownames(cor.mat) = rownames(tmp)
colnames(cor.mat) = rownames(tmp)
cor.mat = melt(cor.mat)
names(cor.mat) = c("protA","protB","rr")
cor.mat$protAB = paste(cor.mat$protA, cor.mat$protB, sep="_")

# n samples in common
nn = nrow(tmp)
n.mat = matrix(rep(NA, nn^2), nrow(tmp), nrow(tmp))
for (ii in 1:nn) {
  print(ii)
  I1 = !is.na(tmp[ii,])
  for (jj in 1:nn) {
    if (ii>=jj) next
    I2 = !is.na(tmp[jj,])
    n.mat[ii,jj] = sum(I1 & I2)
    n.mat[jj,ii] = n.mat[ii,jj]
  }
}
rownames(n.mat) = rownames(tmp)
colnames(n.mat) = rownames(tmp)
n.mat = melt(n.mat)
names(n.mat) = c("protA","protB","nn")
n.mat$protAB = paste(n.mat$protA, n.mat$protB, sep="_")

# TP, TN, or missing
gs.mat = adjacency_matrix_from_list(gold_standard)
a = as.data.frame(t(combn(rownames(tmp), 2)))
protsA = as.character(a[,1])
protsB = as.character(a[,2])
lab_idxs <- protsA %in% rownames(gs.mat) &
  protsB %in% rownames(gs.mat)
idxing_mat <- cbind(protsA[lab_idxs], protsB[lab_idxs])
labels <- rep(NA, length(protsA))
labels[lab_idxs] <- gs.mat[idxing_mat]
labels = data.frame(protAB = paste(protsA, protsB, sep="_"),
                    labels = labels)

df = inner_join(inner_join(cor.mat, n.mat, by="protAB"), labels, by="protAB")

I1 = which(df$nn>4 & df$labels==1)
I0 = which(df$nn>4 & df$labels==0)
df2 = data.frame(rr = c(df$rr[I1], df$rr[I0]), 
                 group = c(rep(1, length(I1)), rep(0, length(I0))))
ggplot(df2, aes(x=rr, fill=factor(group))) + geom_density(alpha=.4)




################################ How good are the complexes
unqprots.data = unique(c(rownames(data[[1]][[1]]), rownames(data[[1]][[2]]), rownames(data[[1]][[3]])))
unqprots.gold = unique(unlist(sapply(gold_standard, strsplit, ";")))

nchan = length(data)
nrep = length(data[[1]])

df.gs = data.frame(prots = character(length(gold_standard)),
                   nn = numeric(length(gold_standard)),
                   nn.overlap = numeric(length(gold_standard)),
                   rrh1 = numeric(length(gold_standard)), 
                   rrh2 = numeric(length(gold_standard)), 
                   rrh3 = numeric(length(gold_standard)), 
                   rrm1 = numeric(length(gold_standard)), 
                   rrm2 = numeric(length(gold_standard)), 
                   rrm3 = numeric(length(gold_standard)), 
                   stringsAsFactors = F)
for (ii in 1:length(gold_standard)) {
  prots = gold_standard[ii]
  
  I.h1 = match(prots, rownames(data[[1]][[1]]))
  
  for (jj in 1:length(prots)) {
    for (kk in 1:length(prots)) {
      
    }
  }
}




