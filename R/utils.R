# function to check if there is an NA in the row
.IsNARow<-function(x){
  if (sum(is.na(x))>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.checkbinary<-function(x){
  if (sum(x==0)+sum(x==1)==length(x) & sum(x==0)>0 & sum(x==1)>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.paircombn<-function(vec1,vec2=NULL,same=TRUE){
  result1<-vector("numeric",0)
  result2<-vector("numeric",0)
  if (same){
    for (i in 1:length(vec1)){
      result1<-c(result1,rep(vec1[i],length(vec1)+1-i))
      result2<-c(result2,vec1[i:length(vec1)])
    }
  } else {
    if (is.null(vec2)){stop("If pairs are not drawn from the same list argument vec2 should be specified.")}
    for (i in 1:length(vec1)){
      result1<-c(result1,rep(vec1[i],length(vec2)))
      result2<-c(result2,vec2)
    }
  }
  lapply(1:length(result1),function(x){c(result1[x],result2[x])})
}












