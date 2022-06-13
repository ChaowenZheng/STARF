
rm(list=ls())
library(optimx,lib.loc = "H:/rlibs/4.1.0")
library(psych)
#library(optimParallel)#for impelmeting parallel estimation
library(ggplot2,lib.loc = "H:/rlibs/4.1.0")
library(ggrepel,lib.loc = "H:/rlibs/4.1.1")
library(Cairo,lib.loc = "H:/rlibs/4.1.0")
# plot of map
library(sf,lib.loc = "H:/rlibs/4.1.0")
library(tmap,lib.loc = "H:/rlibs/4.1.0")
library(sp,lib.loc = "H:/rlibs/4.1.0")
library(tidyverse)

library(maptools)
library(tigris,lib.loc = "H:/rlibs/4.1.0")#for mearging data frame with spatial object
options(tigris_use_cache = TRUE)
library("rnaturalearth",lib.loc = "H:/rlibs/4.1.0")
library("rnaturalearthdata",lib.loc = "H:/rlibs/4.1.0")
library(expm,lib.loc = "H:/rlibs/4.1.0")    ##For calculation of the power of the matrix
library(pheatmap)## for plotting the heat map
#install.packages("farver",dependencies = TRUE)
library(reshape2)##for the melt function

#Setting working directory and put the code "STARF_func_para_first_ando" in this directory
setwd("D:/My Drive/STARF_DATA/Final_tfp/Final code and data/")
source("STARF_func_para_first_ando.R")

# Basic Settings ----------------------------------------------------------

#number of outlier spatial
rmax<-4
Tall<-2

W_range<-c(2)
Wmax<-max(W_range)
Wall<-length(W_range)
cri1<-10^-6
cri2<-9^-7
VIP<-6

#excluding some asia countries and select W>0.01 could make the results better

Exp<-"TFP_Final"
#method<-"BFGS"
method<-"LBFGSB"
Exclude<-"T"
Normalisation<-"T"
# for (set in c("no","two","three")){
for (set in c("one")){
  #Take three values no, two, three
  #number of regressors
  if (set=="one"){
    k<-4
  }
  if (set=="two"){
    k<-5
  }
  if (set=="three"){
    k<-6
  }
  maxite<-50
  summary_sta<-matrix(0,11,k*Wmax*Tall)
  #meaning for each row could be found in the main document
  #Read TFP data
  #dim(tfp)
  #Read data and weighting matrix
  for (tind in 2:Tall){
    if (tind==1){
      data<-read.csv("rgdp_growth_1960_determinants.csv",header=T)
      tfp<-read.csv("tfp_1960.csv")
      dim(tfp)
      tfp<-as.matrix(tfp[,2:ncol(tfp)])
      
    } else {
      data<-read.csv("./data/rgdp_growth_1970_tfp.csv",header=T)
      tfp<-read.csv("./data/tfp_1970.csv")
      tfp<-as.matrix(tfp[,2:ncol(tfp)])
    }
    
    
    for (wind in W_range){
      if (tind==1){
        if (wind==1){
          W<-read.csv("W_inverse_distance2_1960_determinants.csv",header=T)
        }
        if (wind==2){
          W<-read.csv("./data/trade_weights_1960_determinants.csv",header=T)

        }
      } else{
        if (wind==1){
          W<-read.csv("./data/W_inverse_distance2_1970_tfp.csv",header=T)
        }
        if (wind==2){
          W<-read.csv("./data/trade_weights_1970_2015.csv",header=T)

        }
      }
      
      W<-W[2:ncol(W)]
      W<-as.matrix(W)
     
      #Specify the dependent variabe,e.g.

      Y<-as.matrix(data[,4:ncol(data)])
      #make it a N*T matrix
      N<-nrow(Y)
      #number of individuals 2nd row
      summary_sta[2,(tind-1)*Wmax*k+(wind-1)*k+2]<-N
      T<-ncol(Y)
     
      
      #make it a N*T*k matrix
      # In case you want to include a constant term, you need to generate it and then include it in X
      #Why need to indicate the constant term, is because in MC, sometimes, we ignore constant
      #The setting here is more convenient for both empirical analysis and MC simulations
      
      #generate spatial,lagged,spatial lagged variable
      Ylag<-Y[,1:(ncol(Y)-1)]
      Y<-Y[,2:ncol(Y)]
      WY<-W%*%Y
      #is.numeric(Y)
      WYlag<-W%*%Ylag
      dim(WYlag)
      #Specify the exogeneous explanatory variable, including (constant), Y_{t-1} and WY_{t-1}
      #make it a N*T*k matrix
      X_exo<-array(1,dim = c(nrow(Y),ncol(Y), k))
      dim(X_exo)
      X_exo[,,2]<-Ylag
      X_exo[,,3]<-WYlag
      if(set=="one"){
        X_exo[,,4]<-tfp
      }
      if(set=="two"){
        X_exo[,,4]<-rnna
        X_exo[,,5]<-emp
      }
      if(set=="three"){
        X_exo[,,4]<-rnna
        X_exo[,,5]<-emp
        X_exo[,,6]<-tfp
      }
   
      
      
      # Determining the number of factors --------------------------------------------------------------
      
      #Specify the max number of factor, e.g.
      #matrix for storing the information criterion values
      IC_value<-matrix(0,rmax,3) #3 for 3 Informations criterions
      #matrix for storing estimated factors by IC1, IC2 and IC3
      r_hat<-matrix(0,1,3)
      
      for (r_id in 1:rmax){
        #The main working function
        
        ##Input of the function
        #Y: N*T matrix
        #X: N*T*k matrix, where k is the number of regressors and may include constant,Y_{t-1},WY_{t-1}
        #By treating Y_{t-1},WY_{t-1} as regressors, the function is also suitable for cases with other covariates
        #W: N*N spatial weighting matrix
        #R: Number of factor
        #converge_value: stop criterion for stoping iterated estimation
        #max_iter: maximum number of iterations
        #optim_method: the optimization method for numerical serach. Currently, we only provide "Nelder-Mead" and "BFGS" method.
        
        ##Output:
        #theta_hat: the parameters of inetrest
        #F_hat: estimated factors
        #Lambda_hat:estimated factor loadings
        #Varaince-covariance matrix of the estimator.
        #IC_val: information criterion value for determining the number of factors
        
        
        #the following setup is used when running simulation
        # In case the estimation is not satisfactory, you may change them.
        STARF_fit<-starf_function(Y,X_exo,W,r_id,cri1,maxite,method)
        IC_value[r_id,]<-STARF_fit$IC_val
      }
      
      for (ic in 1:3){
        #number of estimated factor,1st row
        summary_sta[1,(tind-1)*Wmax*k+(wind-1)*k+ic]<-which.min(IC_value[,ic])
        r_hat[,ic]<-which.min(IC_value[,ic])
      }
      colnames(r_hat)<-c("IC1","IC2","IC3")
      #Print the estimated unober of factors using the three ICs
      
      
      
      # Estimation --------------------------------------------------------------
      
      #Choose one the selected number of factors,e.g
      #there are 5 optimization mthods that could be used:
      # 1.BFGS 2.Nelder-Mead 3.nlminb 4.newuoa 5.bobyqa
      
      STARF_fit<-starf_function(Y,X_exo,W,r_hat[,2],cri2,maxite,method)
      #STARF_fit<-starf_function(Y,X_exo,W,4,8*10^-6,100,"BFGS")
      
      #Obtain the estimated parameters: N*k matrix, k is the number of interested parameters
      
      theta_hat<-STARF_fit$theta_hat
      
      #Obtain the estimated residuals
      #Notice that res_hat is in general very small since Y itself is very small.
      res_hat<-STARF_fit$residual
      
      
      #Obtain the variance and covariance matrix for the estimators, which is N*k*k array.
      # theta_hat<-matrix(rnorm(N*4),N,4)
      # theta_variance<-matrix(1,N,4)
      theta_variance<-STARF_fit$Variance_hat
      
      theta_variance_final<-matrix(0,N,k+1)
      # Obtain p-value of the estimators
      
      #T statistics
      t_stat<-matrix(0,N,k+1)
      for (i in 1:N){
        theta_variance_final[i,]<-sqrt(diag(theta_variance[i,,])[1:(k+1)])
        t_stat[i,]<-(theta_hat[i,])/theta_variance_final[i,]
      }
      p_value<-(1-pnorm(abs(t_stat)))*2
      # calculation the p values for the estimators: N*K matrix
      sig_num<-matrix(0,N,k+1)
      star<-matrix(0,N,k+1)
      #p_value<-1-pnorm(abs(t_stat))
      for (i in 1:nrow(p_value)){
        for (j in 1:(k+1)){
          if (p_value[i,j]<=0.1){
            sig_num[i,j]<-1
            if(p_value[i,j]<=0.05){
              if(p_value[i,j]<=0.01){
                star[i,j]<-3
              } else{
                star[i,j]<-2
              }
            } else{
              star[i,j]<-1
            }
          } else {
            star[i,j]<-0
          }
        }
      }
      #Number of significant estimates 8 th row
      summary_sta[8,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-colSums(sig_num)[c(1,3:(k+1))]
      # which(abs(theta_hat[,1])>=1)
      # which(abs(p_value[,1])<=0.1)
      #save the estimated results
      results<-cbind(theta_hat,theta_variance_final,p_value,star)
      results<-data.frame(data[,1:3],results)
      write.csv(results,file= sprintf("./results/Individual Estimation Results, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      ##Store the estimated residual for IRF and FEVD analysis
      write.csv(res_hat,file= sprintf("./results/Estimated Residual, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
   
      
      #Read the estimation results
      ####Summary statistics of the estimated results
      theta_hat<-theta_hat[,c(1,3:(k+1))]
      # #number of estimation with abs(spatial)>1,3rd row
      #number of estimation with abs(spatial)>1 hitting the bound,3rd row
      summary_sta[3,((tind-1)*Wmax*k+(wind-1)*k+2)]<-length(which(abs(theta_hat[,1])>=0.996))
      #number of estimation with abs(dynamic)>1,4th row
      summary_sta[4,((tind-1)*Wmax*k+(wind-1)*k+2)]<-length(which(abs(theta_hat[,2])>=1))
      #number of estimation with stationary condition violated, 5th row
      summary_sta[5,((tind-1)*Wmax*k+(wind-1)*k+2)]<-length(which((theta_hat[,1]+theta_hat[,2]+theta_hat[,3])>=1))
      index_unstable<-which((theta_hat[,1]+theta_hat[,2]+theta_hat[,3])>=1)
      #number of positive estimator,6th row
      pos<-matrix(0,1,k)
      #number of negative estimator,7th row
      neg<-matrix(0,1,k)
      for (i in 1:k){
        pos[,i]<-length(which(theta_hat[,i]>0))
        neg[,i]<-length(which(theta_hat[,i]<0))
      }
      summary_sta[6,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-pos
      summary_sta[7,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-neg
      
      
      
      ##Mean group estimator,9th row
      if (Exclude == "T"){
        summary_sta[9,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-colMeans(theta_hat[-index_unstable,])
      } else {
        summary_sta[9,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-colMeans(theta_hat)
      }
      c<-which(rowSums(theta_hat[,1:3])>=1)
      if(length(c)==0){
        summary_sta[10,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-sqrt(diag(var(theta_hat)))
        sd<-sqrt(diag(var(theta_hat)))

        summary_sta[11,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-(1-pnorm(abs(colMeans(theta_hat)/sd)))*2
      } else{
        summary_sta[10,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-sqrt(diag(var(theta_hat[-c,])))
        sd<-sqrt(diag(var(theta_hat[-c,])))
        summary_sta[11,((tind-1)*Wmax*k+(wind-1)*k+1):((tind-1)*Wmax*k+(wind*k))]<-(1-pnorm(abs(colMeans(theta_hat)/sd)))*2
      }
      #Sd of MG estimator, 10th row
      
      write.csv(summary_sta,file= sprintf("./results/Summary/Summary Statistics,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude))
 
      
      ##Density Plot
      # save
      theta_hat<-data.frame(data[,1:3],theta_hat)
      if (Exclude == "T"){
        theta_hat<-theta_hat[-index_unstable,]
      } 
      
      names(theta_hat)[names(theta_hat) == "X1"] <- "Spatial"
      names(theta_hat)[names(theta_hat) == "X2"] <- "Dynamic"
      names(theta_hat)[names(theta_hat) == "X3"] <- "Spatial_Lag"
      names(theta_hat)[names(theta_hat) == "X4"] <- "TFP"
      for (var in c("Spatial","Dynamic","Spatial_Lag","TFP")){
        if(var=="Spatial"){
          density <- ggplot(theta_hat,aes(x=Spatial)) +
            geom_density()+geom_vline(aes(xintercept=mean(Spatial)), color="blue", linetype="dashed", size=1)+
            theme(axis.title.x=element_blank())
        }
        # Basic density
        if(var=="Dynamic"){
          density <- ggplot(theta_hat,aes(x=Dynamic)) +
            geom_density()+geom_vline(aes(xintercept=mean(Dynamic)), color="blue", linetype="dashed", size=1)+
            theme(axis.title.x=element_blank())
        }
        if(var=="Spatial_Lag"){
          density <- ggplot(theta_hat,aes(x=Spatial_Lag)) +
            geom_density()+geom_vline(aes(xintercept=mean(Spatial_Lag)), color="blue", linetype="dashed", size=1)+
            theme(axis.title.x=element_blank())
        }
        if(var=="TFP"){
          density <- ggplot(theta_hat,aes(x=TFP)) +
            geom_density()+geom_vline(aes(xintercept=mean(TFP)), color="blue", linetype="dashed", size=1)+
            theme(axis.title.x=element_blank())
        }
        
        filename <- sprintf('./results/fig/Coefficient_Density_%s,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s',var,N,T,wind,Exp,method,Exclude)
        ggsave(
          filename = sprintf("%s.png", filename),
          plot     = density,
          device   = "png",
          type     = "cairo",
          dpi      = 300)
      }
      
      
      ##Map plot
      world <- ne_countries(scale = "medium", returnclass = "sp")
      #names(world)
      #We need to return sp class, which provides more features
      #we now join the map with our data to selec the countries we use
      world <- geo_join(world, theta_hat, "adm0_a3", "ISO.Code", how = "inner")
      #world$V1
      ##Plot of the estimated coefficients
      tmap_options(check.and.fix = TRUE)
      
      for (var2 in c("Spatial","Dynamic","Spatial_Lag","TFP")){
        
        map<-qtm(world, var2)
        filename <- sprintf('./results/fig/Coefficient_map_%s,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s',var2,N,T,wind,Exp,method,Exclude)
        tmap_save(tm= map, filename = sprintf("%s.png", filename),dpi = 300)
      }
      
      
      #     ###Network Analysis -------------------------------------------------
      
      
      
      
      #     #Individual Analysis ------------------------------------------------
      
      
      ##IRF
      ##The auto-regressive parameter matrix
      N<-nrow(Y)
      if (Exclude=="T"){
        N<-N-length(index_unstable)
        W<-W[-index_unstable,-index_unstable]
        res_hat<-res_hat[-index_unstable,]
      }
      I_N<-diag(1,N,N)
      S<-(I_N-diag(theta_hat$Spatial)%*%W)
      res_ori<-S%*%res_hat
      Psi<-solve(S)%*%(diag(theta_hat$Dynamic)+diag(theta_hat$Spatial_Lag)%*%W)
      IRF<-matrix(0,N,N)
      FEVD_num<-matrix(0,N,N)
      FEVD_de<-matrix(0,N,N)
      FEVD<-matrix(0,N,N)
      DD<-matrix(0,N,N)
      
      ##we should make the correlation to be zero, therefore we only take the diagonal variance
      #Sigma<-var(t(res_ori))
      Sigma<-diag(diag(var(t(res_ori))))
      write.csv(Sigma,file= sprintf("./results/FEVD/Sigma,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      
      #correlation<-cor(t(res_ori))
      horizon<-c(0,1:10)
      #horizon<-c(0)
      FEVD_network<-matrix(0,N,length(horizon)*3)
      DD_network<-matrix(0,N,length(horizon)*3)
      FEVD_vip<-matrix(0,2*N,VIP*length(horizon))
      DD_vip<-matrix(0,2*N,VIP*length(horizon))
      ##TFP
      TFP_effects<-matrix(0,N,N)
      Psi_TFP<-solve(S)%*%(diag(theta_hat$TFP))
      TFP_network<-matrix(0,N,length(horizon)*3)
      TFP_vip<-matrix(0,2*N,VIP*length(horizon))
      
      h_index<-0
      for (h in horizon){
        h_index<-h_index+1
        #we are having an AR(1) process
        Psi_h<-Psi%^%h
        for (i in 1:N){
          e_i<-matrix(0,N,1)
          e_i[i,]<-1
          for (j in 1:N){
            e_j<-matrix(0,N,1)
            e_j[j,]<-1
            # IRF[j,i]<-t(e_j)%*%Psi_h%*%solve(S)%*%Sigma%*%e_i/sqrt((t(e_i)%*%solve(S)%*%Sigma%*%t(solve(S))%*%e_i))
            IRF[j,i]<-t(e_j)%*%Psi_h%*%solve(S)%*%Sigma%*%e_i/sqrt((t(e_i)%*%Sigma%*%e_i))
            FEVD_num[j,i]<-FEVD_num[j,i]+(t(e_j)%*%Psi_h%*%solve(S)%*%Sigma%*%e_i)^2/Sigma[i,i]
            FEVD_de[j,i]<-FEVD_de[j,i]+(t(e_j)%*%(Psi_h)%*%solve(S)%*%Sigma%*%t(solve(S))%*%t(Psi_h)%*%e_j)
            FEVD[j,i]<-FEVD_num[j,i]/FEVD_de[j,i]
          }
        }
        print(h)
        # rowSums(IRF)
        #By construction, the sum of the GFEVD may exceed 1
        ##Row normalization
        # for (i in 1:N){
        #   rowsum<-sum(FEVD[i,])
        #   for (j in 1:N){
        #     FEVD[i,j]<-FEVD[i,j]/rowsum
        #   }
        # }
        #print(rowSums(FEVD))
        write.csv(IRF,file= sprintf("./results/IRF/IRF,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",h,N,T,wind,Exp,method,Exclude),row.names = FALSE)
        write.csv(FEVD,file= sprintf("./results/FEVD/FEVD,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",h,N,T,wind,Exp,method,Exclude),row.names = FALSE)
        #Cumulative Direct effects
        DD<-DD+Psi_h%*%solve(S)
        #write.csv(DD,file= sprintf("./results/Direct_Derivative/DD,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",h,N,T,wind,Exp,method,Exclude),row.names = FALSE)
        
        ##Network Analysis
        
        
        #Cumulative TFP effects
        TFP_effects<-TFP_effects+Psi_h%*%Psi_TFP
        write.csv(TFP_effects,file= sprintf("./results/TFP_effects/TFP_effects,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",h,N,T,wind,Exp,method,Exclude),row.names = FALSE)
        
        for (i in 1:N){
          if (Normalisation=="T"){
            FEVD_network[i,(h_index-1)*3+1]<-sum(FEVD[i,-i])/sum(FEVD[i,])
          } else {
            FEVD_network[i,(h_index-1)*3+1]<-sum(FEVD[i,-i])
          }
          FEVD_network[i,(h_index-1)*3+2]<-sum(FEVD[-i,i])
          FEVD_network[i,(h_index-1)*3+3]<-sum(FEVD[-i,i])-sum(FEVD[i,-i])
          DD_network[i,(h_index-1)*3+1]<-sum(DD[i,-i])
          DD_network[i,(h_index-1)*3+2]<-sum(DD[-i,i])
          DD_network[i,(h_index-1)*3+3]<-sum(DD[-i,i])-sum(DD[i,-i])
          TFP_network[i,(h_index-1)*3+1]<-sum(TFP_effects[i,-i])
          TFP_network[i,(h_index-1)*3+2]<-sum(TFP_effects[-i,i])
          TFP_network[i,(h_index-1)*3+3]<-sum(TFP_effects[-i,i])-sum(TFP_effects[i,-i])
        }
        
        
        
        #Try to extract useful information from the estimated results
        ###Selection of  5 Important countries 
        #5 largest FEVD contributor
        
        for (i in 1:N){
          FEVD_ind<-order(FEVD[i,],decreasing = TRUE)[1:VIP]
          FEVD_vip[(i-1)*2+1,((h_index-1)*VIP+1):(h_index*VIP)]<-data$ISO.Code[FEVD_ind]
          FEVD_vip[(i-1)*2+2,((h_index-1)*VIP+1):(h_index*VIP)]<-round(FEVD[i,FEVD_ind],4)
          DD_ind<-order(abs(DD[i,]),decreasing = TRUE)[1:VIP]
          DD_vip[(i-1)*2+1,((h_index-1)*VIP+1):(h_index*VIP)]<-data$ISO.Code[DD_ind]
          DD_vip[(i-1)*2+2,((h_index-1)*VIP+1):(h_index*VIP)]<-round(DD[i,DD_ind],4)
          TFP_ind<-order(TFP_effects[i,],decreasing = TRUE)[1:VIP]
          TFP_vip[(i-1)*2+1,((h_index-1)*VIP+1):(h_index*VIP)]<-data$ISO.Code[TFP_ind]
          TFP_vip[(i-1)*2+2,((h_index-1)*VIP+1):(h_index*VIP)]<-round(TFP_effects[i,TFP_ind],4)
        }
        
      
        ##Heat map plot
        
        for (var in c("FEVD","DD","TFP")){
          if (var=="FEVD"){
            melted_cormat <- melt(FEVD)
          }
          if (var=="DD"){
            melted_cormat <- melt(DD)
          }
          if (var=="TFP"){
            melted_cormat <- melt(TFP_effects)
          }
          head(melted_cormat)
          
          
          dim(melted_cormat)
          ggheatmap<- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill = value))+
            geom_tile(color = "white")+
            scale_fill_gradientn(values=c(1,0.05,0), colours=c("#006D2C","#E6F5D0","#D1E5F0"),name="Postitive\n\n\n\n\n\n\nZero")+
            
            #colours=c("#006D2C","#4D9221", "#7FBC41", "#B8E186","#E6F5D0","#D1E5F0","#92C5DE","#4393C3","#2171B5", "#08519C"),name="Postitive\n\n\n 0 \n\n\nNegative")+
            #scale_fill_gradient2(low = "blue", high = "green", mid = "white",  midpoint = 0.1,  space = "Lab") +
            theme_minimal()+ 
            theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
            coord_fixed()+ xlab("Countries") + ylab("Countries")+
            theme(
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_text(size=15,face="bold"),
              axis.title.y = element_text(size=15,face="bold"),
              legend.text =element_blank()
            )+guides(fill = guide_colourbar(title.position = "right", title.hjust = .5,label.position = "right"))
          plot(ggheatmap)
          
          filename <- sprintf("./results/fig/heatmap/Heatmap for %s,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",var,h,N,T,wind,Exp,method,Exclude)
          #main<-sprintf("Between Type GCM analysis for %s", variable)
          ggsave(
            filename = sprintf("%s.png", filename),
            plot     = ggheatmap,
            device   = "png",
            type     = "cairo",
            dpi      = 300)
        }
        
        # Regional Anlysis --------------------------------------------------------
        
        regional_infor<-read.csv("./data/country_code_classification.csv",header=T) 
        if(Exclude=="T"){
          regional_infor<-regional_infor[-index_unstable,]
        }
        typelist<-c("Continent","Subregion","Income","Development")
        vlist<-c("FEVD","DD","TFP_effects")
        for (type in typelist ){
          type_element<-unique(regional_infor[,type])
          
          ### Mean Group Estimator
          MG_estimaion<-matrix(0,3*length(type_element),4)
          theta_MG<-theta_hat[,c(5,4,6,7)]
          i<-0
          for (e0 in type_element){
            i<-i+1
            mg_ind<-which(regional_infor[,type]==e0)
            MG_estimaion[(i-1)*3+1,]<-colMeans(theta_MG[mg_ind,])
            if(length(mg_ind)==1){
              MG_estimaion[(i-1)*3+2,]<-0
            } else {
              MG_estimaion[(i-1)*3+2,]<-sqrt(diag(var(theta_MG[mg_ind,])))
            }
            MG_estimaion[(i-1)*3+3,]<-MG_estimaion[(i-1)*3+1,]/MG_estimaion[(i-1)*3+2,]
            MG_estimaion[(i-1)*3+3,]<-(1-pnorm(abs(MG_estimaion[(i-1)*3+3,])))*2
          }
          rownames(MG_estimaion)<-rep(type_element,each=3)
          write.csv(MG_estimaion,file= sprintf("./results/MG/MG Estimation %s, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",type,N,T,wind,Exp,method,Exclude))
          for (v in vlist){
            if(v=="FEVD"){
              effect<-FEVD
            }
            if(v=="DD"){
              effect<-DD
            }
            if(v=="TFP_effects"){
              effect<-TFP_effects
            }
            region_effect<-matrix(0,length(type_element),length(type_element))
            i<-0
            for (e1 in type_element){
              i<-i+1
              i_ind<-which(regional_infor[,type]==e1)
              j<-0
              for(e2 in type_element){
                j<-j+1
                j_ind<-which(regional_infor[,type]==e2)
                region_effect[i,j]<-sum(effect[i_ind,j_ind])/(0.5*(length(j_ind)+length(i_ind)))
              }
            }
            if(v=="FEVD"){
              for (i in length(type_element)){
                fevd_row_sum<-sum(region_effect[i,])
                for(j in length(type_element)){
                  region_effect[i,j]<-region_effect[i,j]/fevd_row_sum
                }
              }
            }
            colnames(region_effect)<-type_element
            rownames(region_effect)<-type_element
            GCM<-matrix(0,5,length(type_element))
            for (i in 1:length(type_element)){
              GCM[1,i]<-sum(region_effect[i,-i])
              GCM[2,i]<-sum(region_effect[-i,i])
              GCM[3,i]<-sum(region_effect[-i,i])-sum(region_effect[i,-i])
            }
            for (i in 1:length(type_element)){
              GCM[4,i]<-GCM[1,i]/sum(abs(region_effect[i,]))
              GCM[5,i]<-GCM[3,i]/(0.5*sum(abs(GCM[3,])))
            }
            rownames(GCM)<-c("RSI","RSO","RNE","EM","SI")
            
            region_effect<-rbind(region_effect,GCM)
            write.csv(region_effect,file= sprintf("./results/GCM/GCM matrix %s,%s,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",type,v,h,N,T,wind,Exp,method,Exclude))
            
            #########2d plot of EM nad SI
            
            em_si<-data.frame(cbind(GCM[4,],GCM[5,]))
            colnames(em_si)<-c("EM","SI")
            em_si$group<-type_element
            
            par(mar=c(1,1,1,1))
            
            
            gcm_plot <- ggplot(em_si, aes(EM, SI, label = group))+xlab("External Motivation") + ylab("Systematic Influence")+
              geom_point(color = ifelse(em_si$SI >=0, "red", "green"))
            
            
            gcm_plot <- gcm_plot + geom_text_repel() #+ labs(title = "geom_text_repel()")
            
            dev.off()
            
            filename <- sprintf("./Results/GCM/GCM_plot/EM_SI plot_%s,%s,h=%d, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",type,v,h,N,T,wind,Exp,method,Exclude)
            ggsave(
              filename = sprintf("%s.png", filename),
              plot     = gcm_plot,
              device   = "png",
              type     = "cairo",
              dpi      = 300)
            # ggsave(
            #   filename = sprintf("%s.pdf", filename),
            #   plot     = gcm_plot,
            #   device   = cairo_pdf)
            }
        }
      }
      
      write.csv(FEVD_vip,file= sprintf("./results/FEVD/FEVD_vip,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      write.csv(FEVD_network,file= sprintf("./results/FEVD/FEVD_network,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      write.csv(TFP_network,file= sprintf("./results/TFP_effects/TFP_network, N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      write.csv(TFP_vip,file= sprintf("./results/TFP_effects/TFP_vip,N=%d,T=%d,W=%d,%s,method=%s,Exclude=%s.csv",N,T,wind,Exp,method,Exclude),row.names = FALSE)
      cat("Experiemnts, T,w, specifciation,Exclude", tind,wind,Exclude, "\n")
    }
  }
}




