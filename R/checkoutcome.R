.checkoutcomes<-function(outcomevector,Force=FALSE,binary=FALSE){
  listoutcomes<-as.numeric(levels(as.factor(outcomevector)))[order(as.numeric(levels(as.factor(outcomevector))))]
  no.outcomes<-length(listoutcomes)
  if (no.outcomes>20 & Force==FALSE){
    stop("More than 20 different values for outcome variable.\n If you are sure you wish to estimate this model rerun command with Force option set to TRUE.")
  }
  if (binary & no.outcomes>2){
    stop("More than 2 values for outcome variable. Try ologit or oprobit as appropriate.")
  }
  output<-list(listoutcomes=listoutcomes,no.outcomes=no.outcomes)
}
