

pancultural <- function(df, group, nfactors) {
  if (nfactors > ncol(df) - 1){
    print("Please specify a lower number of dimensions")
  } else{
  df[,group] <- as.numeric(df[,paste(group)])
  PAN <- psych::statsBy(data = df, group = paste(group),
                        cors = TRUE, cor = "cor", method ="pearson")
  PAN_PCA <- psych::principal(PAN$pooled, nfactors = paste(nfactors),
                              rotate = "varimax")
  return(PAN_PCA$loadings)
  }
}
