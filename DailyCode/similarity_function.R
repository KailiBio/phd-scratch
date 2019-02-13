

####################
# Adjusted Rand Index
####################
ari<-function(x,y)
{
  t=as.matrix(table(x,y));
  TP=sum(choose(t,2));
  FN=sum(choose(apply(t,1,sum),2))-TP;
  FP=sum(choose(apply(t,2,sum),2))-TP;
  TN=sum(choose(sum(t),2))-TP-FN-FP;
  a=TP+1;
  b=TN+1;
  c=FN+1;
  d=FP+1;
  
  ari=(a-(a+c)*(a+d)/(a+b+c+d))/((a+c+a+d)/2-(a+c)*(a+d)/(a+b+c+d));
  return(ari);
}

####################
# Jaccard Index
####################
ji<-function(x,y)
{
  t = unlist(table(apply(cbind(x,y), 1, sum)))
  #
  if(2 %in% names(t)){
    m11 = t[names(t)=="2"]
  }else{
    m11 = 0
  }
  #
  if(1 %in% names(t)){
    m10 = t[names(t)=="1"]
  }else{
    m10 = 0
  }
  if(m10==0 && m11==0){
    ji=-1
  }else{
    ji = m11/(m11+m10)
  }
  return(ji);
}