

source("functions_dani.R")

############################### # plot well-correlated profiles
if (0) {
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
}

################################ # plot OGT interactions
if (1) {
  
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
}



################################ # plot OGA correlators
if (1) {
  oga = "Q9EQQ9"
  ia = which(ml.ints$precision>=0.5 &
               (ml.ints$protein_A==oga | ml.ints$protein_B==oga))
  ib = which(hl.ints$precision>=0.5 &
               (hl.ints$protein_A==oga | hl.ints$protein_B==oga))
  # no oga interactions...
  df = data.frame(fraction = numeric(10^5),
                  value = numeric(10^5),
                  condition  =character(10^5),
                  rep = numeric(10^5), stringsAsFactors = F)
  df2plot = data.frame(fraction = numeric(10^5),
                       value = numeric(10^5),
                       condition  =character(10^5),
                       protein = character(10^5),
                       rep = numeric(10^5), 
                       plot.number = numeric(10^5),stringsAsFactors = F)
  text.layer = data.frame(x = numeric(10^5), y = numeric(10^5),
                          label = character(10^5), 
                          plot.number = numeric(10^5), stringsAsFactors = F)
  df.oga = data.frame(proteinA = character(1000),
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
      ia = which(rownames(tmp)==oga)
      
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
      
      ib = which(cor.mat>.8 & nsamples>quantile(nsamples, 0.8) & nsamples>3)
      for (kk in 1:length(ib)) {
        if (rownames(tmp)[ib[kk]]==oga) next
        cc3 = cc3+1
        I = (cc2+1) : (cc2+60)
        df2plot$fraction[I] = 1:60
        df2plot$value[I] = tmp[ia,]
        df2plot$condition[I] = ii
        df2plot$protein[I] = 1#"Q9EQQ9"
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
        text.layer$label[cc4] = "Q9EQQ9"
        text.layer$plot.number[cc4] = cc3
        
        df.oga$proteinA[cc3] = "Q9EQQ9"
        df.oga$proteinA[cc3] = rownames(tmp)[ib[kk]]
        df.oga$correlation[cc3] = cor.mat[ib[kk]]
        df.oga$condition[cc3] = ii
      }
    }
  }
  df = df[1:cc,]
  df.oga = df.oga[1:cc3,]
  df2plot = df2plot[1:cc2,]
  I = seq(from=1, to=cc4, by=2)
  text.layer = text.layer[1:cc4,]
  text.layer$x = 40
  text.layer$y[I] = 4#max(df.plot$V1, na.rm=T)*.5
  text.layer$y[I+1] = 2.5#max(df.plot$V1, na.rm=T)*.8
  text.layer$protein = ""
  text.layer$protein[I] = 2
  text.layer$protein[I+1] = 1
  text.layer$figure.number = ceiling(text.layer$plot.number/25)
  df2plot$figure.number = ceiling(df2plot$plot.number/25)
  
  for (ii in 1:3) {
    ia = which(df2plot$figure.number == ii)
    ib = which(text.layer$figure.number == ii)
    ggplot(df2plot[ia,], aes(x=fraction, y=value, color=protein)) + facet_wrap(~plot.number) +
      geom_line(alpha=.3) + geom_point(alpha=.5) + 
      geom_text(data=text.layer[ib,], aes(x=x, y=y, label=label), size=3) +
      #geom_text(data=text.layer[ib+1,], aes(x=x, y=y, label=label), size=3) +
      theme_bw() + ylab("protein amount") + theme(strip.text.x = element_blank()) + 
      coord_cartesian(ylim=c(0,5))
    ggsave(paste("../figures/oga_correlators_",ii, ".png", sep=""), width=8, height=4.8)
  }
  
  # just plot oga
  df$condition[df$condition==1] = "M/L"
  df$condition[df$condition==2] = "H/L"
  df$rep = as.character(df$rep)
  ggplot(df, aes(x=fraction,y=value,color=rep,group=rep)) + 
    geom_point(alpha=.5) +geom_line(alpha=.5) +
    facet_wrap(~condition)
  ggsave("../figures/oga_bycondition.png", width=6, height=2.8)
}




################################ # plot all interactions
# this one is big!!
if (0) {
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
  for (ii in 1:length(I2plot)) {
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
    sf = paste("/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/figures/interactions/", 
               ii, "_", protA, "_", protB, ".png", sep="")
    ggsave(sf, width=8, height=3.5, dpi = 100)
  }
}



################################ # mike's functional analysis
if (1) {
  # read ontology
  ontology = get_ontology("../data/go-basic.obo")
  # read annotations
  goa = read_gpa("../data/mgi.gpad.gz")
  # remove roots 
  rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
  goa %<>% dplyr::filter(!GO.ID %in% rootNames)
  
  # map MGI to Uniprot
  mgi2uni = as.data.frame(read_tsv("../data/mgi_to_uniprot.txt"))
  ia = match(goa$UNIPROT, mgi2uni$`yourlist:M201911268471C63D39733769F8E060B506551E1258D67FC`)
  goa$UNIPROT = as.character(mgi2uni$Entry[ia])
  
  # process BP, CC, and MF annotations separately
  bp = filter_roots(goa, ontology, 'BP') %>% 
    as_annotation_list("UNIPROT", "GO.ID")
  cc = filter_roots(goa, ontology, 'CC') %>% 
    as_annotation_list("UNIPROT", "GO.ID")
  mf = filter_roots(goa, ontology, 'MF') %>% 
    as_annotation_list("UNIPROT", "GO.ID")
  anns = list(BP = bp, CC = cc, MF = mf)
  # create overall annotation object
  #ann = as_annotation_list(goa, "UNIPROT", "GO.ID")
  
  df = differential2(ml.ints, hl.ints, anns, ontology, nboot = 1000)
  
  # check all interactions are in interactome
  for (ii in 1:nrow(df)) {
    tmpml = unlist(strsplit(df$interactions.ml[ii], ";"))
    tmphl = unlist(strsplit(df$interactions.hl[ii], ";"))
    if (length(tmpml)>0){
      for (jj in 1:length(tmpml)) {
        prots = sort(unlist(strsplit(tmpml[jj], "-")))
        if (!sum(ml.ints$protein_A ==prots[1] & ml.ints$protein_B==prots[2])==1) error
      }
    }
    if (length(tmphl)>0) {
      for (jj in 1:length(tmphl)) {
        prots = sort(unlist(strsplit(tmphl[jj], "-")))
        if (!sum(hl.ints$protein_A ==prots[1] & hl.ints$protein_B==prots[2])==1) error
      }
    }
  }
  
  # write
  write_tsv(df, path = "../data/dani-differential.txt")
}