print.margins.oglmx<-function(x, ... ){
  for (m in 1:length(x)){
    cat("Marginal Effects on Pr(Outcome==",names(x)[m],")","\n",sep="")
    if (m==length(x)){
      printCoefmat(x[[m]])
    } else {
      printCoefmat(x[[m]],signif.legend=FALSE)
      cat("------------------------------------","\n")
    }
  }
}
