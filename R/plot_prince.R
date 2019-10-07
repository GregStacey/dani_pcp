
# plot chromatograms
interactions.hl$unique.interaction = character(nrow(interactions.hl))
for (ii in 1:nrow(interactions.hl)) {
  interactions.hl$unique.interaction[ii] = paste(sort(interactions.hl[ii,1:2]), collapse="_")
}

interactions.ml$unique.interaction = character(nrow(interactions.ml))
for (ii in 1:nrow(interactions.ml)) {
  interactions.ml$unique.interaction[ii] = paste(sort(interactions.ml[ii,1:2]), collapse="_")
}

interactions = full_join(interactions.hl, interactions.ml, by="unique.interaction")

# interactions in both
xx = 0.5
I = which(interactions$precision.x>xx & interactions$precision.y>xx)

# interactions in just ML
I1 = which(interactions$precision.x>xx & interactions$precision.y<xx)
I2 = which(interactions$precision.x<xx & interactions$precision.y>xx)
labels = c("M+H interaction", "Only M interaction", "Only H interaction")

I2plot = c(I, I1, I2)
for (ii in 1:length(I2plot)) {
  protA = interactions$protein_A.x[I2plot[ii]]
  protB = interactions$protein_B.x[I2plot[ii]]
  
  ia = which(rownames(dani[[1]]) == protA)
  ib = which(rownames(dani[[1]]) == protB)
  
  if (I2plot[ii] %in% I) {
    label = labels[1]
  } else if (I2plot[ii] %in% I1){
    label = labels[2]
  } else if (I2plot[ii] %in% I2){
    label = labels[3]
  }
  
  df = data.frame(protA.ml = dani[[1]][ia,],
                  protB.ml = dani[[1]][ib,],
                  protA.hl = dani[[2]][ia,],
                  protB.hl = dani[[2]][ib,], stringsAsFactors = F)
  df$fraction = as.numeric(sapply(sapply(as.character(rownames(df)), strsplit, ".", fixed=T), "[", 2))
  df = melt(df, id.vars = "fraction")
  df$variable = as.character(df$variable)
  df$channel = sapply(sapply(df$variable, strsplit, ".", fixed=T), "[", 2)
  df$prot = sapply(sapply(df$variable, strsplit, ".", fixed=T), "[", 1)
  df$prot[df$prot=="protA"] = protA
  df$prot[df$prot=="protB"] = protB
  ggplot(df, aes(x=fraction, y=value, color=prot)) + geom_line() + geom_point() +
    facet_wrap(~channel) + ggtitle(label)
  sf = paste("/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/figures/chromatograms/", 
             ii, "_", protA, "_", protB, ".png", sep="")
  ggsave(sf, width=10, height=3.5, dpi = 100)
}
