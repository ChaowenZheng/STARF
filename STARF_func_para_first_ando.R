
# Intruduction ------------------------------------------------------------

# This Programme is for estimating a heterogeneous spatial dynamic model with unobserved factors

# The algorithms:
####Step 1: Obtain Initial estimaion for all the parameters
#1). Uisng OLS to obtain the initial estimator for parameters of interest by ignoring the factors
#2) Apply PC/SVD to the residual sequence ofr initial estimation of factor/loadings

####Step 2: Update the estimator
#3) Based on all the other estimators, update the spatial coefficient according to (6) in Ando and Lu(2019) using numerical search
#4) Based on the updated spatial coefficient, initial factor/loadings, updated the other parameters of interest by OLS
#5) Update factor/loadings by PC/SVD to the residual series based on the updated parehters of  interest
#6) Repeat Step 3-5 until convergence.


# Estimation function -------------------------------------------------------------

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


starf_function<-function(Y,X,W,R,converge_value, max_iter,optim_method){
  
  # Baisc Setting -------------------------------------------------------
  # X<-X_exo
  # R<-2
  # converge_value<-5*10^-7
  # max_iter<-100
  # optim_method<-"LBFGSB"
  N<-nrow(Y)
  I_N<-diag(1,N,N)
  T<-ncol(Y)
  k1<-dim(X)[3]
  #Number of all explanatory variable (plus the spatial lagged variable)
  k2<-k1+1
  #Generate spatial veriable
  WY<-W%*%Y
  
  #STAGE 1: Generate intial value ---------------------------------------------------
  
  
  #Step1: Initial value of parameters by OLS, named theta_ini
  theta_ini<-matrix(0,N,k2)
  resid_ini<-matrix(0,N,T)
  for (i in 1:N) {
    X_i<-cbind(WY[i,],X[i,,])
    fit<-lm(Y[i,]~X_i-1)
    theta_ini[i,]<-fit$coefficients
    resid_ini[i,]<-fit$residuals
  } 
  
  
  #Step2: Apply SVD to the residual matrix to obtain initial factor and factor loadings
  #Normalization are satisfied automatically
  F_ini<- sqrt(T)*svd(resid_ini)$v[,1:R]
  #The following is the same
  # VEC <- eigen(t(resid_ini)%*%resid_ini)$vectors; 
  # F_ini2 <- sqrt(T)*(VEC)[,1:R]; 
  Lambda_ini<-resid_ini%*%F_ini%*%solve(t(F_ini)%*%F_ini)
  LF_ini<-Lambda_ini%*%t(F_ini)
  
  mean_ini<-colMeans(theta_ini)
  cat("Initial Mean of Parameters= ", mean_ini, "\n")
  # #Better Initial
  #
  # # It is obvious that the above initial value is biased
  # # we then iterate the above procedure to obtain a better initial value
  # # This is mainly for speeding up the computation:
  # # When updating the parameters, it would take more time for the result to converge starting from a bad initial value.
  # # By contrast, iterated OLS estimation is much quicker.
  # # Of course, this step could be removed if you have any other concern.
  
  
  diff<-1 #Differnence between two estimates
  theta_update<-matrix(0,N,k2)
  U_ini<-matrix(0,N,T)#U=LF+epsilon
  Niter<-0
  while (diff > 10^-6&& Niter< 100) {
    
    #Step 1': potential better initial value for parameters
    for (i in 1:N) {
      X_i<-cbind(WY[i,],X[i,,])
      fit<-lm((Y[i,]-LF_ini[i,])~X_i-1)
      theta_update[i,]<-fit$coefficients
      U_ini[i,]<-Y[i,]-X_i%*%theta_update[i,]
    }
    
    ##Step 2': potential better initial value for factor/loadings
    ##apply SVD to y-rhoWY-X'beta,Normalization are satisfied automatically
    F_update<- sqrt(T)*svd(U_ini)$v[,1:R]
    # VEC <- eigen(t(U_ini)%*%U_ini)$vectors;
    # F_update <- sqrt(T)*(VEC)[,1:R];
    Lambda_update<-U_ini%*%F_update %*%solve(t(F_update)%*%F_update)
    LF_update<-Lambda_update%*%t(F_update)
    
    #Comparision of the new and old estimated values
    diff<-mean((theta_update-theta_ini)^2)+mean((LF_update-LF_ini)^2)
    Niter<-Niter+1
    
    theta_ini<-theta_update
    F_ini<-F_update
    Lambda_ini<-Lambda_update
    LF_ini<-LF_update
    # mean_ini<-colMeans(theta_ini)
    # cat("Initial Mean of Parameters= ", mean_ini,diff, "\n")
  }
  mean_ini<-colMeans(theta_ini)
  cat("Initial Mean of Parameters= ", mean_ini, "\n")
  
  
  # STAGE 2: Updating the estimators -------------------------------------------------
  
  
  diff<-1 #Differnence between two estimates
  theta_update<-matrix(0,N,k2)
  Niter<-0
  while (diff > converge_value && Niter< max_iter) {
    
    #Step3:Update spatial Coefficient
    #by using the numerical serach to find the minmizer of  (Y-(I-rhoW)^{-1}(LF+X'beta))^2
    #based on initial values of factor/loading/ and other parameters of interest.
    
    X_update_spatial<-LF_ini
    for(p in 1:k1){
      X_update_spatial<-X_update_spatial+diag(theta_ini[,1+p])%*%X[,,p]
    }
    lb<-rep(-0.996,N)
    ub<-rep(0.996,N)
    ###Using numerical serach to find the minimizer since the objective function is nonlinaer in parameters of interest
    #The input is N*1 vector of estimated heterogeneous spatial coefficient.
    if(optim_method=="P_BFGS"||optim_method=="P_LBFGSB" ){
      if(optim_method=="P_BFGS"){
      objective_func<-function(rho,Y,W,X_update_spatial,N,T) {
        #we restrict the spatial coefficient to be less than 1 for more precise estimation
        rho[which(abs(rho)>=1)]<-0.996
        Mrho = diag(rho)
        N<-nrow(Mrho)
        I_N<-diag(1,N,N)
        obj_value<-psych::tr((Y-solve(I_N-Mrho%*%W)%*%X_update_spatial)%*%t(Y-solve(I_N-Mrho%*%W)%*%X_update_spatial))/(N*T)
        return(obj_value)
      }
      } else {
        objective_func<-function(rho,Y,W,X_update_spatial,N,T) {
          #we restrict the spatial coefficient to be less than 1 for more precise estimation
          #rho[which(abs(rho)>1)]<-0.996
          Mrho = diag(rho)
          N<-nrow(Mrho)
          I_N<-diag(1,N,N)
          obj_value<-psych::tr((Y-solve(I_N-Mrho%*%W)%*%X_update_spatial)%*%t(Y-solve(I_N-Mrho%*%W)%*%X_update_spatial))/(N*T)
          return(obj_value)
        }
      }
    } else { 
      if(optim_method=="LBFGSB"){
      objective_func<-function(rho) {
        #we restrict the spatial coefficient to be less than 1 for more precise estimation
        Mrho = diag(rho)
        N<-nrow(Mrho)
        I_N<-diag(1,N,N)
        obj_value<-psych::tr((Y-solve(I_N-Mrho%*%W)%*%X_update_spatial)%*%t(Y-solve(I_N-Mrho%*%W)%*%X_update_spatial))/(N*T)
        return(obj_value)
      }} else {
        objective_func<-function(rho) {
          #we restrict the spatial coefficient to be less than 1 for more precise estimation
          rho[which(abs(rho)>=1)]<-0.996
          Mrho = diag(rho)
          N<-nrow(Mrho)
          I_N<-diag(1,N,N)
          obj_value<-psych::tr((Y-solve(I_N-Mrho%*%W)%*%X_update_spatial)%*%t(Y-solve(I_N-Mrho%*%W)%*%X_update_spatial))/(N*T)
          return(obj_value)
        }
      }
    }
    #There are many numerical search method.and we provide the following two options for this.
    #For more information, see  help(optimx)
    if(optim_method=="Nelder-Mead"){
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("Nelder-Mead"), control = list(reltol = 1e-7, abstol = 1e-7))
    }
    if(optim_method=="BFGS"){
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("BFGS"), control = list(reltol = 1e-7, abstol = 1e-7))
    }
    if(optim_method=="LBFGSB"){
      #It is very important that the intial value should satisfy the constraints
      theta_ini[which(abs(theta_ini[,1])>0.996),1]<-0.996
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("L-BFGS-B"), control = list(reltol = 1e-7, abstol = 1e-7),lower=lb,upper=ub)
    }
    if(optim_method=="P_BFGS"){
      #use the parallel optimization procedure to save time
      spatial_fit <- optimParallel(par=theta_ini[,1],Y=Y,W=W,X_update_spatial=X_update_spatial,N=N,T=T,fn=objective_func,  method = "BFGS")
    }
    if(optim_method=="P_LBFGSB"){
      #use the parallel optimization procedure to save time
      spatial_fit <- optimParallel(par=theta_ini[,1],Y=Y,W=W,X_update_spatial=X_update_spatial,N=N,T=T,fn=objective_func,  method = "L-BFGS-B",lower=lb,upper=ub)
    }
    if(optim_method=="nlminb"){
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("nlminb"), control = list(reltol = 1e-7, abstol = 1e-7))
    }
    if(optim_method=="newuoa"){
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("newuoa"), control = list(reltol = 1e-7, abstol = 1e-7))
    }
    if(optim_method=="bobyqa"){
      spatial_fit= optimx(theta_ini[,1], objective_func, method = c("bobyqa"), control = list(reltol = 1e-7, abstol = 1e-7))
    }
    #relative tolerance and absolute tolrence
    if(optim_method=="P_BFGS"||optim_method=="P_LBFGSB"){
      theta_update[,1]<-as.matrix(spatial_fit$par)
    }
    else{
      theta_update[,1]<-as.matrix(spatial_fit[1:N])
    }
    
    #Step 4: Update parameter for exogeneous regressors, including constant,Y_{t-1},WY_{t-1}
    #Via using the OLS regression for y-rhoWY-LF on X
    
    Y_update_exo<-Y-LF_ini-diag(theta_update[,1])%*%WY
    resid_update<-matrix(0,N,T)
    for (i in 1:N) {
      X_i<-X[i,,]
      fit<-lm(Y_update_exo[i,]~X_i-1)
      theta_update[i,2:k2]<-fit$coefficients
    } 
    
    #Step 5: Update factors/loadings
    resid_update<-Y-diag(theta_update[,1])%*%WY
    for(p in 1:k1){
      resid_update<-resid_update-diag(theta_update[,1+p])%*%X[,,p]
    }
    ##apply SVD to y-rhoWY-X'beta,Normalization are satisfied automatically
    F_update<- sqrt(T)*svd(resid_update)$v[,1:R]
    # VEC <- eigen(t(resid_update)%*%resid_update)$vectors; 
    # F_update <- sqrt(T)*(VEC)[,1:R];
    Lambda_update<-resid_update%*%F_update %*%solve(t(F_update)%*%F_update)
    LF_update<-Lambda_update%*%t(F_update)
    
    #Comparision of the intial and updated estimation
    diff<-mean((theta_update-theta_ini)^2)+mean((LF_update-LF_ini)^2)
    
    
    #Number of iterations
    Niter<-Niter+1
    Mean_final<-colMeans(theta_ini)#-colMeans(theta0)

    theta_ini<-theta_update
    F_ini<-F_update
    Lambda_ini<-Lambda_update
    LF_ini<-LF_update
    cat("Iteration, Mean of Estimators, Diff = ", Niter,Mean_final,diff, "\n")
  }
  theta_final<-theta_update
  F_final<-F_update
  Lambda_final<-Lambda_update
  
  
  # Variance calculation ----------------------------------------------------
  cat("Claculating Variance...","\n")
  
  #Collect all the regrssores as X_star
  psi_star<-array(0,dim = c(N,T, 1+k1+R))
  #Calculate the transformed spatial variable
  X_sum<-Lambda_final%*%t(F_final)
  for(p in 1:k1){
    X_sum<-X_sum+diag(theta_final[,1+p])%*%X[,,p]
  }
  #Calculate the matrix (I-rhoW)^{-1}
  W_star<-solve(I_N-diag(theta_final[,1])%*%W)
  WY_star<-W%*%W_star%*%X_sum
  psi_star[,,1]<-WY_star
  
  #Calculation of the variance of the idiosyncratic error
  e_hat<-Y-W_star%*%X_sum
  e_hat2<-rowMeans(e_hat^2)
  
  #Collect the other variables
  for(p in 1:k1){
    psi_star[,,1+p]<-X[,,p]
  }
  F<-array(0,dim = c(N,T,R))
  for(i in 1:N){
    F[i,,]<-F_final
  }
  for(p in 1:R){
    psi_star[,,1+k1+p]<-F[,,p]
  }
  
  #Matrix for storing variance
  variance_all<-array(0,dim=c(N,k1+1+R,k1+1+R))
  #variance_all2<-array(0,dim=c(N,k1+1+R,k1+1+R))
  for(k in 1:N){
    #Calculation of Gamma_k and Phi_k
    Gamma_k<-matrix(0,1+k1+R,1+k1+R)
    Phi_k<-matrix(0,1+k1+R,1+k1+R)
    #n_i<-matrix(0,N,1)
    for(i in 1:N){
      # for (t in 1:T){
      #   Gamma_k<-Gamma_k+W_star[i,k]^2*as.matrix(psi_star[k,t,])%*%psi_star[k,t,]
      #   Phi_k<-Phi_k+e_hat2[i]*W_star[i,k]^2*as.matrix(psi_star[k,t,])%*%(psi_star[k,t,])
      # }
      Gamma_k<-Gamma_k+W_star[i,k]^2*t(psi_star[k,,])%*%psi_star[k,,]
      Phi_k<-Phi_k+e_hat2[i]*W_star[i,k]^2*t(psi_star[k,,])%*%(psi_star[k,,])
    }
    variance_all[k,,]<-solve(Gamma_k)%*%Phi_k%*%solve(Gamma_k)
  }
  
  ##Information Criterion value for determining the number of factors
  #See equation (9) of Bai and Ng(2002)
  #Following Bai and Li (2021), the penalty function is divided by 2
  #Notice that since Ando and Lu(2019) use a reduced model for estimaton, the error term is also transformed.
  ICvalue<-matrix(0,1,3)
  
  ICvalue[,1]=log(mean(e_hat^2))+R*(N+T)*log(N*T/(N+T))/(N*T)
  ICvalue[,2]=log(mean(e_hat^2))+R*(N+T)*log(min(N,T))/(N*T)
  ICvalue[,3]=log(mean(e_hat^2))+R*log(min(N,T))/(min(N,T))
  
  #Store all the estimation results
  list(theta_hat=theta_final,F_hat=F_final,Lambda_hat=Lambda_final,IC_val=ICvalue,Variance_hat=variance_all,residual=e_hat)
}




