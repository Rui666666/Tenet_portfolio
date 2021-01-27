rm(list = ls(all = TRUE))
graphics.off()

wdir = "/Users/LvB/Documents/GitHub/Tenet-portfolio"
setwd(wdir)

library(quadprog)
library(lpSolve)
library(stringr)
library(ggplot2)
library(dplyr)
library(Matrix)

source("Portfolio_methods.R")

#----------------------------------------Data preparation ----------------------------------------

# data period
date_start_data = 20180701 
date_end_data = 20200924  

# porfolio construction 
date_start = 20190106
date_end = 20200924


#price
price = read.csv("Crypto_Price_20200924.csv", header = TRUE)
lambda = read.csv("lambdas_wide_20180701_20200924_90_0.05.csv", header = TRUE)
tau=0.05

price[is.na(price)] = 0
price$date=as.numeric(gsub("-", "", price$date))

return=diff(log(as.matrix(price[,-1]))) 
date=price$date[-1]
return[is.na(return)]=0
for (i in (1:nrow(return))){
  for (j in (1:ncol(return))){
    if ((return[i,j]=="Inf")|| (return[i,j]=="-Inf")){
      return[i,j]=0
    }
  }
}

return=cbind(date,return)
return=as.data.frame(return)
return=return[return$date>=date_start_data,]
price=price[price$date>=date_start_data,]
lambda=lambda[lambda$date>=date_start_data,]

lambda$date==return$date
lambdas=lambda[,c(1:17)]  
return=return[colnames(lambdas)]
price=price[colnames(lambdas)]


#----------------------------------------Application----------------------------------------

## 1. Comulatative wealth
ticker = return$date
N = nrow(return)       # number of total days
N0 = which( ticker == date_start )
N1 = which( ticker == date_end )
s = 90                #rolling window lengthn / days in the sample / construction period
reb = 1              #rebalancing period
N_upd = floor((N1-N0)/reb+1)  
Er =0.008             #expected return  mean(mu)
nss = 0               #Equals 1 if no short sales, for markowitz
ga=0.8               #qtec gamma 

pweight_m = matrix(0,N_upd,ncol(return)) #markowitz
wealth_m = matrix(1,N_upd,2)
pweight_q = matrix(0,N_upd,ncol(return)) #qtec
wealth_q = matrix(1,N_upd,2)
pweight_l = matrix(0,N_upd,ncol(return)) #ltec
wealth_l = matrix(1,N_upd,2)


## for ltec
for (t in seq(N0, N1, reb)){
  rdata = return[(t-s+1):t, ]
  ldata = lambdas[(t-s+1):t, ]
  #squre root
  ldata[,2:ncol(lambdas)] = sqrt(lambdas[(t-s+1):t, 2:ncol(lambdas)] ) 
  
  #portfolio parameter
  mu =colMeans(rdata[,c(2:ncol(rdata))])    #mean vector
  n = length(mu)                            #number of assests
  cv = cov(rdata[,-1])                      #return covariance matrix
  l_cv=cov(ldata[,-1])                      #lambda covariance matrix
  l_mu=colMeans(ldata[,c(2:ncol(ldata) )])  #lambda mean
  sol_l = ltec(n,mu,l_mu,Er)             
  w_l = as.matrix(sol_l$solution)
  pweight_l[(t-N0)/reb+1,1] = ticker[t]
  pweight_l[(t-N0)/reb+1,c(2:ncol(return))] = w_l
  #portfolio wealth
  if (t==N0){
    wealth_l[(t-N0)/reb+1,1] = ticker[t]    
  } else {
    ret = log(price[t,-1]) - log(price[t-reb,-1])   
    for (i in (1:length(ret))){
      if ( (ret[i]=="-Inf") || (ret[i]=="NaN") || (ret[i]=="Inf")){
        ret[i]=0
      }
    }
    wealth_l[(t-N0)/reb+1,1] = ticker[t]           
    wealth_l[(t-N0)/reb+1,2] = wealth_l[(t-N0)/reb,2] + as.numeric(ret) %*% pweight_l[(t-N0)/reb,-1]
  }
}

#qtec
wealth_q = matrix(1,N_upd,2)
for (t in seq(N0, N1, reb)){
  rdata = return[(t-s+1):t, ]
  ldata = lambdas[(t-s+1):t, ]  
  #portfolio parameter
  mu =colMeans(rdata[,c(2:ncol(rdata))])    
  n = length(mu)                   
  cv = cov(rdata[,-1])             
  l_cv=cov(ldata[,-1])             
  l_mu=colMeans(ldata[,c(2:ncol(ldata) )])   
  
  w_q =qtec(n,mu,l_mu,l_cv,Er,nss,ga)
  pweight_q[(t-N0)/reb+1,1] = ticker[t]
  pweight_q[(t-N0)/reb+1,c(2:ncol(return))] = w_q
  if (t==N0){
    wealth_q[(t-N0)/reb+1,1] = ticker[t]    
  } else {
    ret = log(price[t,-1]) - log(price[t-reb,-1])   
    for (i in (1:length(ret))){
      if ( (ret[i]=="-Inf") || (ret[i]=="NaN") || (ret[i]=="Inf")){
        ret[i]=0
      }
    }
    wealth_q[(t-N0)/reb+1,1] = ticker[t]           
    wealth_q[(t-N0)/reb+1,2] = wealth_q[(t-N0)/reb,2] + as.numeric(ret) %*% pweight_q[(t-N0)/reb,-1]
  }
}

## markowitz
for (t in seq(N0, N1, reb)){
  rdata = return[(t-s+1):t, ]
  ldata = lambdas[(t-s+1):t, ]
  #portfolio parameter
  mu =colMeans(rdata[,c(2:ncol(rdata))])    
  n = length(mu)                   
  cv = cov(rdata[,-1])            
  l_cv=cov(ldata[,-1])            
  l_mu=colMeans(ldata[,c(2:ncol(ldata) )])   
  # portfolio weight
  sol_m = markowitz_nss(n,mu,cv,Er,nss)  
  w_m = as.matrix(sol_m$solution)        
  pweight_m[(t-N0)/reb+1,1] = ticker[t]
  pweight_m[(t-N0)/reb+1,c(2:ncol(return))] = w_m
  #portfolio wealth
  if (t==N0){
    wealth_m[(t-N0)/reb+1,1] = ticker[t]     
  } else {
    ret = log(price[t,-1]) - log(price[t-reb,-1])   
    for (i in (1:length(ret))){
      if ( (ret[i]=="-Inf") || (ret[i]=="NaN") || (ret[i]=="Inf")){
        ret[i]=0
      }
    }
    wealth_m[(t-N0)/reb+1,1] = ticker[t]           
    wealth_m[(t-N0)/reb+1,2] = wealth_m[(t-N0)/reb,2] + as.numeric(ret) %*% pweight_m[(t-N0)/reb,-1]    #Wt=W(t-1)+w(t-1)Xt
  }
}


# #----------------------------------------Portfolio result----------------------------------------
# plot wealth for different strategies
df<- data.frame(
  x=as.Date(as.character(wealth_l[,1]),  "%Y%m%d"),
  y1=wealth_m [,2],
  y2=wealth_l [,2],
  y3=wealth_q [,2]
)

# Plot
p=ggplot(df, aes(x=x)) +
  scale_x_date(date_labels = "%Y%m%d",breaks='1 month',expand = c(0, 0)) +
  geom_line(aes(y = y1), color="dimgrey") +
  geom_line(aes(y = y2), color="steelblue") +
  geom_line(aes(y = y3), color="tomato2") +
  theme( panel.grid=element_blank(),
         panel.border = element_rect(color="black",fill = NA),
         legend.position = "none",
         axis.title = element_blank(),
         panel.background = element_rect(fill = "transparent",colour = NA),
         plot.background = element_rect(fill = "transparent",colour = NA),
         legend.box.background = element_rect(fill = "transparent"),
         #axis.line = element_line(colour = "black"),
         axis.text.x = element_text(angle = 90, hjust = 1,size=14, face= "bold"),
         axis.text.y = element_text(size=14, face= "bold" ) )

png(paste0("fig_Wealth_",date_end,"_",s, "_",reb, "_",Er,"_",ga,"_",tau,".png"),width = 800, height = 500,
    bg = "transparent")  
p
dev.off()

