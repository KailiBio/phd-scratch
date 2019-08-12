
setwd("/Users/Kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/data_figures/state_jaccard/")

###################
# function
###################
correlation<-function(x,y)
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
  mcc=(a*b-d*c)/sqrt((a+d)*(a+c)*(b+d)*(b+c));
  
  t2 = unlist(table(apply(cbind(x,y), 1, sum)))
  #
  if(2 %in% names(t2)){
    m11 = t2[names(t2)=="2"]
  }else{
    m11 = 0
  }
  #
  if(1 %in% names(t2)){
    m10 = t2[names(t2)=="1"]
  }else{
    m10 = 0
  }
  if(m10==0 && m11==0){
    ji=-1
  }else{
    ji = m11/(m11+m10)
  }
  #jac= a/(a+c+d);
  return(c(ari,mcc,ji));
}

######################
# 1. for same list
######################
pdf("perfect_score_for_JI_ARI_MCC.pdf")
plot(NA, xlim=c(0,100), ylim=c(0,1), xlab="size", ylab = "JI/ARI/MCC")
legend("bottomright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"),
       bty="n")
for(i in 2:100){
  m = as.integer(i*0.5)
  x = c(rep(1,m),rep(0,(i-m)))
  y = c(rep(1,m),rep(0,(i-m)))
  r = correlation(x,y)
  points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
  points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
  points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
}
dev.off()


######################
# 2. ARI & MCC
######################

pdf("similarity_of_ARI_MCC.pdf")
plot(NA, xlim=c(-1,1), ylim=c(-1,1), xlab="ARI", ylab="MCC")
lines(c(0,1), c(0,1), lwd=2, lty=2, col="grey")
for(i in 1:10000){
  for(j in 1:10){
    x = sample(c(0,1), replace=TRUE, size=i)
    y = sample(c(0,1), replace=TRUE, size=i)
    r = correlation(x,y)
    points(r[1],r[2], col = rgb(255,0,0, max=255,alpha=100), pch=20)
  }
}
dev.off()


######################
# 3. differ size, 50% 1
######################

pdf("differ_size_50percent_3methods.pdf")
plot(NA, xlim=c(10,1000), ylim=c(-1,1), xlab="size", ylab="JI & ARI & MCC",
     main = "1: 50% & 0: 50%\nramdom order")
legend("topright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"))
for(i in seq(10,10000,25)){
  for(j in 1:10){
    x = sample(0:1, i, replace=T,prob=c(0.5,0.5))
    y = sample(0:1, i, replace=T,prob=c(0.5,0.5))
    r = correlation(x,y)
    #
    points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
    points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
    points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
  }
}
dev.off()


pdf("differ_size_10percent_3methods.pdf")
plot(NA, xlim=c(10,1000), ylim=c(-1,1), xlab="size", ylab="JI & ARI & MCC",
     main = "1: 10% & 0: 90%, ramdom order")
legend("topright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"))
for(i in seq(10,10000,25)){
  for(j in 1:10){
    m = as.integer(0.1*i)
    n = i-m
    #
    x = sample(0:1, i, replace=T,prob=c(0.9,0.1))
    y = sample(0:1, i, replace=T,prob=c(0.9,0.1))
    r = correlation(x,y)
    #
    points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
    points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
    points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
  }
}
dev.off()

pdf("differ_size_1percent_3methods.pdf")
plot(NA, xlim=c(100,1000), ylim=c(-1,1), xlab="size", ylab="JI & ARI & MCC",
     main = "1: 1% & 0: 99%, ramdom order")
legend("topright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"))
for(i in seq(100,10000,25)){
  for(j in 1:10){
    x = sample(0:1, i, replace=T,prob=c(0.99,0.01))
    y = sample(0:1, i, replace=T,prob=c(0.99,0.01))
    r = correlation(x,y)
    #
    points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
    points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
    points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
  }
}
dev.off()

######################
# 4. length=1000, change state_proportion
######################
pdf("percentage_of_positive_JI_ARI_MCC.pdf")
plot(NA, xlim=c(0,1), ylim=c(-1,1), xlab="percentage of positive states(1)", ylab="JI & ARI & MCC",
     main = "similarity score with percentage of positive values\nsize=1000")
legend("bottomright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"))
for(i in seq(0,1,0.01)){
  for(j in 1:10){
    x = sample(0:1, 1000, replace=T,prob=c(1-i,i))
    y = sample(0:1, 1000, replace=T,prob=c(1-i,i))
    r = correlation(x,y)
    #
    points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
    points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
    points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
  }
}
dev.off()

######################
# 5. perfect score for ARI&MCC in different state_proportion
######################
pdf("percentage_of_positive_perfect_score_for_JI_ARI_MCC.pdf")
plot(NA, xlim=c(0,1), ylim=c(-1,1), xlab="percentage of positive states(1)", ylab="JI & ARI & MCC",
     main = "similarity score in different percentage of positive values\nidentical list,size=1000")
legend("bottomright",legend = c("JI","ARI","MCC"), pch = c(20,17,18), col = c("blue","red","green"),
       bty="n")
for(i in seq(0,1,0.01)){
  for(j in 1:10){
    x = sample(0:1, 1000, replace=T,prob=c(1-i,i))
    y = x
    r = correlation(x,y)
    #
    points(i,r[1], col = rgb(255,0,0, max=255,alpha=100), pch=17)
    points(i,r[2], col = rgb(0,255,0, max=255,alpha=100), pch=18)
    points(i,r[3], col = rgb(0,0,255, max=255,alpha=100), pch=20)
  }
}
dev.off()


######################
