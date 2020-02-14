
require(reshape2)
require(ggplot2)
require(readr)
require(PrInCE)
require(dplyr)
require(tidyr)
require(magrittr)
require(progress)
require(flavin)
require(ontologyIndex)
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
fn = "../data/condition1.csv"
tmp = as.data.frame(read_csv(fn))
tmp = tmp[!grepl("__", tmp$protein.ids), ] # remove REV and CON
#tmp$protein.ids = sapply(sapply(tmp$protein.ids, strsplit, "|", T), "[", 2)
condition1 = (as.matrix(tmp[,3:ncol(tmp)]))
rownames(condition1) = tmp$protein.ids
condition1 = condition1[!is.na(rownames(condition1)),]
condition1 = list(condition1[tmp$replicate==1,], condition1[tmp$replicate==2,], condition1[tmp$replicate==3,])
fn = "../data/condition2.csv"
tmp = as.data.frame(read_csv(fn))
tmp = tmp[!grepl("__", tmp$protein.ids), ] # remove REV and CON
#tmp$protein.ids = sapply(sapply(tmp$protein.ids, strsplit, "|", T), "[", 2)
condition2 = (as.matrix(tmp[,3:ncol(tmp)]))
rownames(condition2) = tmp$protein.ids
condition2 = condition2[!is.na(rownames(condition2)),]
condition2 = list(condition2[tmp$replicate==1,], condition2[tmp$replicate==2,], condition2[tmp$replicate==3,])
data = list(condition1, condition2)

# load gold_standard
fn = "../data/allComplexesmapped.txt"
tmp = as.data.frame(read_tsv(fn))
gs.mapped = sapply(tmp$subunits.UniProt.IDs., strsplit, ";")
names(gs.mapped) = tmp$ComplexName
fn = "../data/allComplexes.txt"
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



# write good interaction lists
fn = "../data/interactions/interactions_reps7_gs1.Rda"
load(fn)
sf = "../data/interactions/interactions_ML.txt"
write_tsv(ml.ints[which(ml.ints$precision>=0.5),], path=sf)
sf = "../data/interactions/interactions_HL.txt"
write_tsv(hl.ints[which(hl.ints$precision>=0.5),], path=sf)


################################ Diagnostics

################################ Are TP and TN pairs different?

#how many proteins in common b/w gold stanard and data?
unqprots.data = unique(c(rownames(data[[1]][[1]]), rownames(data[[1]][[2]]), rownames(data[[1]][[3]])))
unqprots.gold = unique(unlist(sapply(gold_standard[[1]], strsplit, ";")))
n.overlap = length(intersect(unqprots.gold, unqprots.data))
print(paste(length(unqprots.data),"unique proteins in data"))
print(paste(length(unqprots.gold),"unique proteins in gold standard"))
print(paste(n.overlap, " proteins in common"))

# how correlated are true positives? true negatives? random?
uu = 1 # condition
mm = 2 # replicate

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
gs.mat = adjacency_matrix_from_list(gold_standard[[2]])
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




################################ plot well-correlated profiles
conds = list("ml", "hl")
for (uu in 1:2) {
  for (mm in 1:3) {
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
    
    I = which(n.mat$nn>25 & cor.mat$rr>.9 & cor.mat$rr<1)
    
    for (ii in 1:20) {
      ia = which(rownames(data[[uu]][[mm]]) == (n.mat$protA[I[ii]]))
      ib = which(rownames(data[[uu]][[mm]]) == (n.mat$protB[I[ii]]))
      df = data.frame(fraction = rep(1:60, 2),
                      ratio = c(data[[uu]][[mm]][ia,], data[[uu]][[mm]][ib,]),
                      protein = c(rep(as.character(n.mat$protA[I[ii]]), 60), 
                                  rep(as.character(n.mat$protB[I[ii]]), 60)),
                      stringsAsFactors = F)
      ggplot(df, aes(x=fraction, y=ratio, color=protein)) + geom_line() + geom_point()
      sf = paste("../figures/chromatograms/cond",conds[[uu]],"_rep",
                 mm,"_",n.mat$protA[I[ii]],"_",n.mat$protB[I[ii]],".png", sep="")
      ggsave(sf, width=6, height=3)
    }
  }
}



################################ How good are the complexes
gs = gold_standard[[2]] # mapped
unqprots.data = unique(c(rownames(data[[1]][[1]]), rownames(data[[1]][[2]]), rownames(data[[1]][[3]])))
unqprots.gold = unique(unlist(sapply(gs, strsplit, ";")))

nchan = length(data)
nrep = length(data[[1]])

df.gs = data.frame(prots = character(length(gs)),
                   nn = numeric(length(gs)),
                   nn.overlap = numeric(length(gs)),
                   rrh1 = rep(NA,length(gs)), 
                   rrh2 = rep(NA,length(gs)), 
                   rrh3 = rep(NA,length(gs)), 
                   rrm1 = rep(NA,length(gs)), 
                   rrm2 = rep(NA,length(gs)), 
                   rrm3 = rep(NA,length(gs)), 
                   stringsAsFactors = F)
for (ii in 1:length(gs)) {
  prots = unlist(gs[ii])
  df.gs$prots[ii] = paste(prots, collapse=";")
  
  df.gs$nn[ii] = length(prots)
  
  I.m1 = match(prots, rownames(data[[1]][[1]])) %>% na.omit
  I.m2 = match(prots, rownames(data[[1]][[2]])) %>% na.omit
  I.m3 = match(prots, rownames(data[[1]][[3]])) %>% na.omit
  I.h1 = match(prots, rownames(data[[2]][[1]])) %>% na.omit
  I.h2 = match(prots, rownames(data[[2]][[2]])) %>% na.omit
  I.h3 = match(prots, rownames(data[[2]][[3]])) %>% na.omit

  if (length(I.m1)>=2) df.gs$rrm1[ii] = as.data.frame(as.vector(cor(t(data[[1]][[1]][I.m1,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
  if (length(I.m2)>=2) df.gs$rrm2[ii] = as.data.frame(as.vector(cor(t(data[[1]][[2]][I.m2,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
  if (length(I.m3)>=2) df.gs$rrm3[ii] = as.data.frame(as.vector(cor(t(data[[1]][[3]][I.m3,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
  if (length(I.h1)>=2) df.gs$rrh1[ii] = as.data.frame(as.vector(cor(t(data[[2]][[1]][I.h1,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
  if (length(I.h2)>=2) df.gs$rrh2[ii] = as.data.frame(as.vector(cor(t(data[[2]][[2]][I.h2,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
  if (length(I.h3)>=2) df.gs$rrh3[ii] = as.data.frame(as.vector(cor(t(data[[2]][[3]][I.h3,]), , use = "p"))) %>% 
    rename(x = 1) %>% filter(!x==1) %>% unlist %>% mean
}

I = df.gs$nn>3
cor(df.gs[I,4:9], use = "p", method = 'spearman')


################################ 
fns = dir("../data/interactions/", full.names = T)
nn = 15000 
df = data.frame(rank = numeric(nn* 14 * 2),
                score = numeric(nn* 14 * 2),
                precision = numeric(nn* 14 * 2),
                channel = character(nn* 14 * 2),
                protA = numeric(nn* 14 * 2),
                protB = numeric(nn* 14 * 2),
                gs = numeric(nn* 14 * 2),
                rep = numeric(nn* 14 * 2), stringsAsFactors = F)
cc = 0
for (ii in 1:length(fns)) {
  load(fns[ii])
  
  gs = as.numeric(gsub("gs","", gsub(".Rda","", unlist(strsplit(basename(fns[ii]), "_"))[3])))
  rep = as.numeric(gsub("reps","", gsub(".Rda","", unlist(strsplit(basename(fns[ii]), "_"))[2])))
  rep = paste(reps[[rep]], collapse=",")
  
  tmp = list(ml.ints, hl.ints)
  conds = c("ml", "hl")
  for (jj in 1:length(tmp)) {
    I = (cc+1) : (cc+nn)
    df$rank[I] = 1:nn
    df$score[I] = tmp[[jj]]$score[1:nn]
    df$precision[I] = tmp[[jj]]$precision[1:nn]
    df$channel[I] = conds[jj]
    df$protA[I] = tmp[[jj]]$protein_A[1:nn]
    df$protB[I] = tmp[[jj]]$protein_B[1:nn]
    df$gs[I] = gs
    df$rep[I] = rep
    cc = cc+nn
  }
  
  ia = which(ml.ints$precision>0.5 | is.na(ml.ints$precision))
  ib = which(hl.ints$precision>0.5 | is.na(hl.ints$precision))
  print(paste(ii, ": M/L interactions=", 
              length(ia),
              ",  H/L interactions=",
              length(ib), 
              ",  total interactions=",length(ia)+length(ib), sep=""))
}

ggplot(df, aes(x=rank, y=precision, color = channel)) + 
  facet_grid(gs ~ rep) + 
  geom_line()
ggsave("../figures/interaction_scan.png", width = 12, height = 3)



################################ 
# do functional differential analysis on best-performing interactions
# reps = c(1,2,3), gs = un-mapped
# THIS IS NOW DEPRECATED
#   needs to be made compatible with differential-explore.R

#source("/differential.R")
#source("/Users/gregstacey/Academics/Foster/PCP-SILAC/differential/R/differential.R")

# get interactions
fn = "../data/interactions/interactions_reps7_gs1.Rda"
load(fn)

# # filter to 50% precision
# I = which(hl.ints$precision>=0.5)
# nn = I[length(I)]
# hl.ints = hl.ints[1:nn, ]
# I = which(ml.ints$precision>=0.5)
# nn = I[length(I)]
# ml.ints = ml.ints[1:nn, ]
# 
# # fill in precision=NA with precision=1
# I = which(!is.na(hl.ints$precision))[1]
# hl.ints$precision[1:(I-1)] = 1
# I = which(!is.na(ml.ints$precision))[1]
# ml.ints$precision[1:(I-1)] = 1

# ERROR
# For some reason, the first few thousand interactions look bad.
# Therefore, only use interactions with !is.na() precision.
ml.ints = ml.ints[ml.ints$precision>=0.5 & !is.na(ml.ints$precision), ]
hl.ints = hl.ints[hl.ints$precision>=0.5 & !is.na(hl.ints$precision), ]

# run functional differential
dfsig = differential(ml.ints, hl.ints)

# write functional differential analysis
nn = sum(!is.na(dfsig$difference.zscore))
I = order(dfsig$difference.zscore, decreasing = T)
fn = "../data/functional_differential.txt"
write_tsv(dfsig[I[1:nn], ], path=fn)

# attach interactions, proteins, and GO IDs to "../data/functional_differential.txt"
unqprots = "S4R2S6"
for (ii in 1:2) {
  for (jj in 1:3) {
    unqprots = unique(c(unqprots,  rownames(data[[ii]][[jj]])))
  }
}
dfsig = as.data.frame(read_tsv("../data/dani-differential.txt"))
ia = match(dfsig$go.name, ontology$name)
dfsig$go.id = ontology$id[ia]
dfsig$proteins = sapply(dfsig$go.id, FUN=function(x) {
  paste(unlist(ann[[x]][ann[[x]] %in% unqprots]), collapse = " ")
})
dfsig$interactions.hl = character(nrow(dfsig))
dfsig$interactions.ml = character(nrow(dfsig))
for (ii in 1:nrow(dfsig)) {
  prots = unlist(strsplit(dfsig$proteins[ii], " "))
  I.ml = which(ml.ints$protein_A %in% prots & ml.ints$protein_B %in% prots & !is.na(ml.ints$precision))
  I.hl = which(hl.ints$protein_A %in% prots & hl.ints$protein_B %in% prots & !is.na(ml.ints$precision))
  dfsig$interactions.ml[ii] = paste(paste(ml.ints$protein_A[I.ml], 
                                          ml.ints$protein_B[I.ml], sep="-"), collapse = " ; ")
  dfsig$interactions.hl[ii] = paste(paste(hl.ints$protein_A[I.hl], 
                                          hl.ints$protein_B[I.hl], sep="-"), collapse = " ; ")
}
write_tsv(dfsig, path="../data/dani-differential-2.txt")



################################ 
# plot OGT interactions
ogt = "Q8CGY8"
ia = which(ml.ints$precision>=0.5 &
             (ml.ints$protein_A==ogt | ml.ints$protein_B==ogt))
ib = which(hl.ints$precision>=0.5 &
             (hl.ints$protein_A==ogt | hl.ints$protein_B==ogt))
# no OGT interactions...
df = data.frame(fraction = numeric(10^4),
                value = numeric(10^4),
                condition  =character(10^4),
                rep = numeric(10^4), stringsAsFactors = F)
df2plot = data.frame(fraction = numeric(10^4),
                     value = numeric(10^4),
                     condition  =character(10^4),
                     protein = character(10^4),
                     rep = numeric(10^4), 
                     plot.number = numeric(10^4),stringsAsFactors = F)
text.layer = data.frame(x = numeric(10^4), y = numeric(10^4),
                        label = character(10^4), 
                        plot.number = numeric(10^4), stringsAsFactors = F)
df.ogt = data.frame(proteinA = character(1000),
                    proteinB = character(1000),
                    correlation = numeric(1000), 
                    condition = character(1000),stringsAsFactors = F)
cc = 0
cc2 = 0
cc3 = 0
cc4 = 0
for (ii in 1:2) { # conditions
  for (jj in 1:3) { # reps
    tmp = data[[ii]][[jj]]
    ia = which(rownames(tmp)==ogt)
    
    I = (cc+1) : (cc+60)
    df$fraction[I] = 1:60
    df$value[I] = tmp[ia,]
    df$condition[I] = ii
    df$rep[I] = jj
    cc = cc+60

    cor.mat = cor(t(tmp), use = "p")
    cor.mat = cor.mat[ia,]

    # n samples
    nsamples = numeric(nrow(tmp))
    for (kk in 1:nrow(tmp)) {
      nsamples[kk] = sum(!is.na(tmp[ia,]) & !is.na(tmp[kk,]))
    }
    
    ib = which(cor.mat>.8 & nsamples>10)
    for (kk in 1:length(ib)) {
      if (rownames(tmp)[ib[kk]]==ogt) next
      cc3 = cc3+1
      I = (cc2+1) : (cc2+60)
      df2plot$fraction[I] = 1:60
      df2plot$value[I] = tmp[ia,]
      df2plot$condition[I] = ii
      df2plot$protein[I] = 1#"Q8CGY8"
      df2plot$rep[I] = jj
      df2plot$plot.number[I] = cc3
      cc2 = cc2+60
      I = (cc2+1) : (cc2+60)
      df2plot$fraction[I] = 1:60
      df2plot$value[I] = tmp[ib[kk],]
      df2plot$condition[I] = ii
      df2plot$protein[I] = 2#rownames(tmp)[ib[kk]]
      df2plot$rep[I] = jj
      df2plot$plot.number[I] = cc3
      cc2 = cc2+60
      
      cc4 = cc4+1
      text.layer$label[cc4] = rownames(tmp)[ib[kk]]
      text.layer$plot.number[cc4] = cc3
      cc4 = cc4+1
      text.layer$label[cc4] = "Q8CGY8"
      text.layer$plot.number[cc4] = cc3
      
      df.ogt$proteinA[cc3] = "Q8CGY8"
      df.ogt$proteinA[cc3] = rownames(tmp)[ib[kk]]
      df.ogt$correlation[cc3] = cor.mat[ib[kk]]
      df.ogt$condition[cc3] = ii
    }
  }
}
df = df[1:cc,]
df.ogt = df.ogt[1:cc3,]
df2plot = df2plot[1:cc2,]
I = seq(from=1, to=cc4, by=2)
text.layer = text.layer[1:cc4,]
text.layer$x = 40
text.layer$y[I] = 4#max(df.plot$V1, na.rm=T)*.5
text.layer$y[I+1] = 2.5#max(df.plot$V1, na.rm=T)*.8
text.layer$protein = ""
text.layer$protein[I] = 1
text.layer$protein[I+1] = 2

ggplot(df2plot, aes(x=fraction, y=value, color=protein)) + facet_wrap(~plot.number) +
  geom_line(alpha=.3) + geom_point(alpha=.5) + 
  geom_text(data=text.layer[I,], aes(x=x, y=y, label=label), size=3, color="#00BFC4") +
  geom_text(data=text.layer[I+1,], aes(x=x, y=y, label=label), size=3, color="red") +
  theme_bw() + ylab("protein amount") +theme(strip.text.x = element_blank()) + 
  coord_cartesian(ylim=c(0,5))
ggsave("../figures/ogt_correlators.png", width=8, height=4.8)

# just plot ogt
df$condition[df$condition==1] = "M/L"
df$condition[df$condition==2] = "H/L"
df$rep = as.character(df$rep)
ggplot(df, aes(x=fraction,y=value,color=rep,group=rep)) + 
  geom_point(alpha=.5) +geom_line(alpha=.5) +
  facet_wrap(~condition)
ggsave("../figures/ogt_bycondition.png", width=6, height=2.8)

