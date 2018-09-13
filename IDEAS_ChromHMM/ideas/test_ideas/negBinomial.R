
# -- Kaili
# This script is for normalizing HM data by negative binomial distribution (yu's function).

normdata<-function(case, control, method="nbinom")
{
  d=dnorm(-10:10,0,10);
  u= c(rep(0,10),control,rep(0,10));
  l=length(case);
  nu=rep(0,l);
  for(i in 1:21){
    nu = nu + u[i-1+1:l] * d[i];
  }
  nu = nu / sum(d);

  t = which(case<quantile(case,prob=0.99))
  p=mean(case[t])/var(case[t]);
  if(p>0.9) p=0.9;
  if(p<0.1) p=0.1;
  n=mean(case[t])*p/(1-p);
  print(c(n,p));
  #nc = log2((case*20+1)/(nu/mean(nu)*mean(case)*20+1)+1)-1;
  if(method=="nbinom"){
    nc = -pnbinom(case,size=c((nu+1)/mean(nu+1)*n),prob=p,lower.tail=F,log.p=T)/log(10);
  } 
  else{
    nc = -ppois(case,(nu+1)/mean(nu+1)*mean(case),lower.tail=F,log.p=T)/log(10);
    }
  return(nc);
}


args <- commandArgs(trailingOnly=TRUE);
caseFile = args[1]
controlFile = args[2]
outFile = args[3]

case = read.table(caseFile, header = FALSE)
control = read.table(controlFile, header = FALSE)

out = normdata(case[,5], control[,5], method="nbinom")

write.table(out, outFile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
