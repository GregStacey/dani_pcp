
# read proteingroups
#fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/1streplicate_proteinGroups_Dani.txt"
fns = c("../data/1stR_proteinGroups.txt","../data/2ndR_proteinGroups.txt","../data/3rdR_proteinGroups.txt")

condition1 = list() # M/L
condition2 = list() # H/L
for (uu in 1:length(fns)) {
  fn = fns[uu]
  tmp = as.data.frame(read_tsv(fn))
  
  names(tmp) = gsub("-","_", names(tmp))
  
  # extract columns
  # H
  I.hl = which(grepl("Ratio H/L normalized ", names(tmp)))
  # order columns by fraction
  fraction = as.numeric(gsub("Ratio H/L normalized ", "", sapply(strsplit(names(tmp)[I.hl], "_"), "[", 1)))
  I.hl = I.hl[ order(fraction, decreasing = F)]
  # add protein IDs
  I.hl = c(1, I.hl)
  
  # M
  I.ml = which(grepl("Ratio M/L normalized ", names(tmp)))
  # order columns by fraction
  fraction = as.numeric(gsub("Ratio M/L normalized ", "", sapply(strsplit(names(tmp)[I.ml], "_"), "[", 1)))
  I.ml = I.ml[ order(fraction, decreasing = F)]
  # add protein IDs
  I.ml = c(1, I.ml)
  
  # make condition1
  condition1[[uu]] = tmp[,I.hl]
  # make replicate column
  condition1[[uu]]$replicate = uu
  condition1[[uu]] = condition1[[uu]][,c(1,62,2:61)]
  names(condition1[[uu]]) = c("protein.ids", "replicate", paste("fraction", 1:60))
  
  # make condition2
  condition2[[uu]] = tmp[,I.ml]
  # make replicate column
  condition2[[uu]]$replicate = uu
  condition2[[uu]] = condition2[[uu]][,c(1,62,2:61)]
  names(condition2[[uu]]) = c("protein.ids", "replicate", paste("fraction", 1:60))
}
condition1 = rbind(condition1[[1]], condition1[[2]], condition1[[3]])
condition2 = rbind(condition2[[1]], condition2[[2]], condition2[[3]])

# write 
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/condition1"
write_tsv(condition1, paste(fn, ".txt", sep=""))
write_csv(condition1, paste(fn, ".csv", sep=""))
fn = "/Users/gregstacey/Academics/Foster/LabMembers/Dani/pcp/data/condition2"
write_tsv(condition2, paste(fn, ".txt", sep=""))
write_csv(condition2, paste(fn, ".csv", sep=""))


