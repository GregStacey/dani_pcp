
source("functions_dani.R") # ml.ints, hl.ints
ml.ints = ml.ints[ml.ints$precision >= 0.5 & !is.na(ml.ints$precision), ]
hl.ints = hl.ints[hl.ints$precision >= 0.5 & !is.na(hl.ints$precision), ]

# plot chromatograms
hl.ints$unique.interaction = character(nrow(hl.ints))
for (ii in 1:nrow(hl.ints)) {
  hl.ints$unique.interaction[ii] = paste(sort(hl.ints[ii,1:2]), collapse="_")
}

ml.ints$unique.interaction = character(nrow(ml.ints))
for (ii in 1:nrow(ml.ints)) {
  ml.ints$unique.interaction[ii] = paste(sort(ml.ints[ii,1:2]), collapse="_")
}

# order by max precision
interactions = full_join(hl.ints, ml.ints, by="unique.interaction")
tmp = cbind(interactions$precision.x,interactions$precision.y)
I = order(apply(tmp, 1, max, na.rm=T), decreasing = T)
interactions = interactions[I,]

# interactions in both
xx = 0.5
I = which(interactions$precision.x>xx & interactions$precision.y>xx)

# interactions in just ML
I1 = which(interactions$precision.x>xx & is.na(interactions$precision.y))
I2 = which(interactions$precision.y>xx & is.na(interactions$precision.x))
labels = c("M+H interaction", "Only M interaction", "Only H interaction")

I2plot = c(I, I1, I2)
conds = c("M/L", "H/L")
for (ii in 8548:length(I2plot)) {
  prots = unlist(strsplit(interactions$unique.interaction[I2plot[ii]], "_"))
  protA = prots[1]
  protB = prots[2]
  
  if (I2plot[ii] %in% I) {
    label = labels[1]
  } else if (I2plot[ii] %in% I1){
    label = labels[2]
  } else if (I2plot[ii] %in% I2){
    label = labels[3]
  }
  
  df = data.frame(protA = numeric(1e4), 
                  protB = numeric(1e4),
                  channel = numeric(1e4),
                  replicate = numeric(1e4),
                  fraction = numeric(1e4),stringsAsFactors = F)
  cc = 0
  for (jj in 1:length(data)) { #channel
    for (kk in 1:length(data[[jj]])) { #replicate
      ia = which(rownames(data[[jj]][[kk]]) == protA)
      ib = which(rownames(data[[jj]][[kk]]) == protB)
      i0 = (cc+1) : (cc+60)
      df$fraction[i0] = 1:60
      df$protA[i0] = data[[jj]][[kk]][ia,]
      df$protB[i0] = data[[jj]][[kk]][ib,]
      df$channel[i0] = conds[jj]
      df$replicate[i0] = kk
      cc = cc+length(i0)
    }
  }
  df = df[1:cc,]
  df = melt(df, measure.vars = c("protA", "protB"), 
            id.vars = c("protA", "protB", "replicate", "channel", "fraction"))
  df$variable = as.character(df$variable)
  df$variable[df$variable=="protA"] = protA
  df$variable[df$variable=="protB"] = protB
  
  ggplot(df, aes(x=fraction, y=value, color=variable)) + geom_line() + geom_point() +
    facet_grid(channel~replicate) + ggtitle(label) + theme_bw()
  sf = paste("/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/figures/chromatograms/3reps/", 
             ii, "_", protA, "_", protB, ".png", sep="")
  ggsave(sf, width=8, height=3.5, dpi = 100)
}
