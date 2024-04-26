
#####################################################
##### Newton-Raphson iteration: h-quasi-likelihood with lognormal frailty ######

HQL_lnkm=function(beta_h0ln,v_h0ln,alpha_h0ln){
  
  ########################################################
  ##### PART 1: adjust h-quasi-likelihood function: y(alpha) to estimate the frailty component alpha ########
  ###### lognormal frailty ######
  halphaln=function(alpha){
    
    i_alpln=1/alpha;
    v_sqln=sum(v_h1ln*v_h1ln);  ## (v_hat)^2
    ##### 2nd part in h: log density of u ########
    temp2ln=-0.5*(n*log(alpha)+i_alpln*v_sqln);
    
    
    ##### 1st part in h: log density of x ########
    temp1_0ln=colSums(M_01km*cv_1ln);  ### D*1 vector
    temp1_1ln=colSums(M_01km)*c(log(co_Wi1ln[1:(D-num_pr0km)]),rep(0,num_pr0km));  ### D*1 vector
    temp1ln=sum((temp1_0ln-temp1_1ln)*d_n);
    
    
    ######## Hessian matrix of beta ##########
    ####### initial vectors #######
    Hv_xs0ln=NA;Hbv_xsln=NA;
    
    for (a in 1:n) {
      if(a>1){
        n_0=sum(q[1:(a-1)])+1; ### the start number of cluster a
        n_1=sum(q[1:a]);   ### the end number of cluster a
      }else{
        n_0=1;n_1=q[1];
      }
      xi_v2ln=colSums(xi_v21ln*dMc_hat[,n_0:n_1]);  
      xi_xv3ln=-v_h1ln[a]*i_alpln;
      ######### final estimator of xi_v ###########
      xi_xvln=xi_v1ln[n_0:n_1]+xi_v2ln+xi_xv3ln;    
      ######## estimator of Sigma_v[a] ############
      Hv_xs0ln[a]=sum(xi_xvln*xi_xvln);   
      ######## estimator of Sigma_bv[a] ############
      Hbv_xsln[a]=sum(xi_bln[n_0:n_1]*xi_xvln);   
    }
    
    ######## final estimator of Sigma_v: n*n matrix ############
    ### log-normal frailty
    Hv_xsln=diag(Hv_xs0ln);    
    
    ############# the part of -d^2l2/dv^2: n*n matrix ##############
    ### log-normal frailty
    h_xv31ln=diag(i_alpln,n,n); 
    ############# final estimate of A_v: -d^2l/dv^2 ##############
    Hv_xa1ln=h_v11ln+h_v21ln+h_xv31ln;  
    
    ### built Sigma matrix ###
    H_xs1ln=c(Hb_sln,Hbv_xsln);  ### p*n matrix
    H_xs2ln=cbind(Hbv_xsln,Hv_xsln);  ### n*p matrix
    H_xsln=rbind(H_xs1ln,H_xs2ln);  ### (p+n)*(p+n) matrix
    
    
    ##################### update inverse Hessian matrix ######################
    ###############################     
    ### log-normal frailty
    hiv_xa1ln=solve(Hv_xa1ln); ### the inverse matrix of A_v
    H1_xco1ln=colSums(Hbv_a1ln*hiv_xa1ln);  ## n*1 vector
    H2_xco1ln=Hb_a1ln-sum(Hbv_a1ln*H1_xco1ln);
    #### the 1st part in inverse Hessian matrix ##
    H11_x1ln=1/H2_xco1ln; 
    #### the 2nd part in inverse Hessian matrix ##
    H12_x1ln=-H1_xco1ln*H11_x1ln; ### n*1 vector
    #### the 4th part in inverse Hessian matrix ##
    H22_x1ln=hiv_xa1ln+H1_xco1ln%*%t(H1_xco1ln)/H2_xco1ln; ###  n*n matrix
    
    
    ### built inverse Hessian matrix ###
    H_x1ln=c(H11_x1ln,H12_x1ln);  ### p*n matrix
    H_x2ln=cbind(H12_x1ln,H22_x1ln); ### n*p matrix
    H_xln=rbind(H_x1ln,H_x2ln);  ### (p+n)*(p+n) matrix
    
    ### covariance matrix for theta=(beta,v) ###
    eiv_ln=eigen(H_xln)$values;
    
    ######## the adjustment term ##########
    temp3ln=0.5*(p+n)*log(2*pi)+0.5*sum(log(eiv_ln));  
    
    ##### final result #####
    templn=temp1ln+temp2ln+temp3ln;
    return(-templn);
  }
  
  
  
  ######### iteration starts here  ##########
  iter=1;div=0;   ### iteration and diverging times
  repeat
  {  
    ########### log-normal frailty  #########
    v_0ln=rep(v_h0ln,q); ### N*1 vector
    cv_0ln=x1*beta_h0ln+v_0ln;    ### eta_hat: N*1 vector
    expeta_h0ln=as.numeric(exp(cv_0ln));  ### exp(eta_hat): N*1 vector
    ########### common parts include beta and v #########
    Wi_exp0ln=Wi_1km*expeta_h0ln;  #### N*D matrix
    Wx_exp0ln=Wi_x1km*expeta_h0ln;  #### N*D matrix
    co_Wx0ln=colSums(Wx_exp0ln);  #### D*1 vector
    co_Wi0ln=colSums(Wi_exp0ln);  #### D*1 vector
    co_Wxs0ln=colSums(Wi_x1skm*expeta_h0ln);  #### D*1 vector
    
    
    ##### initial matrix: numerator of v_bar #####
    v_mat0ln=matrix(0,nrow=D,ncol=n);  
    vz_mat0ln=matrix(0,nrow=D,ncol=n);  
    for (a in 1:n) {
      if(a>1){
        n_0=sum(q[1:(a-1)])+1; ### the start number of cluster a
        n_1=sum(q[1:a]);   ### the end number of cluster a
      }else{
        n_0=1;n_1=q[1];
      }
      v_mat0ln[,a]=colSums(Wi_exp0ln[n_0:n_1,]);  #### D*1 vector
      vz_mat0ln[,a]=colSums(Wx_exp0ln[n_0:n_1,]);  #### D*1 vector  
    }
    
    #### estimate m0, v_bar and x_bar #######  
    if(num_pr0km>0){
      ## since sum(co_Wx0[,wk_0km])=0, we just consider (D-wk_0km) vector
      mh_wk0=rep(0,num_pr0km);  ## wk_0km*1 vector
      x_wk0=rep(x_mean,num_pr0km); ## wk_0km*1 vector
      v_wk0=matrix(1/n,nrow=num_pr0km,ncol=n); ## wk_0km*n matrix 
      
      ###### log-normal frailty
      mh_0ln=(Mc_0km/co_Wi0ln)[-wk_0km]; ## (D-wk_0km)*1 vector 
      v_bar0ln=v_mat0ln[-wk_0km,]/co_Wi0ln[-wk_0km]; ## (D-wk_0km)*n matrix 
      x_bar0ln=(co_Wx0ln/co_Wi0ln)[-wk_0km]; ## (D-wk_0km)*1 vector
      ### final estimator of m0, v_bar and x_bar: let m0[t(D)]=0,v_bar[t(D)]=1/n,x_bar[t(D)]=x_mean
      mh_0ln=c(mh_0ln,mh_wk0); #### D*1 vector 
      v_bar0ln=rbind(v_bar0ln,v_wk0); #### D*n matrix 
      x_bar0ln=c(x_bar0ln,x_wk0); #### D*1 vector
    }else{ 
      
      ### estimator of m0[t(k)] and x_bar[t(k)]: k=1,...,D ############
      mh_0ln=Mc_0km/co_Wi0ln; ## D*1 vector 
      v_bar0ln=v_mat0ln/co_Wi0ln; ## D*n matrix
      x_bar0ln=co_Wx0ln/co_Wi0ln; ## D*1 vector
    }
    
    
    #### estimate M(t(k)): N*D matrix #####
    M_02ln=t(t(Wi_exp0ln)*mh_0ln);
    M_hat0ln=M_01km-M_02ln;
    
    
    ############a) estimating equation of beta #########
    sc_bln=colSums(M_hat0ln*x1); #### D*1 vector
    dft_bln=sum(sc_bln*d_n); 
    
    
    ########### Hessian matrix of beta ############
    ######## estimate of A_b ############
    Hb_a0ln=sum((co_Wxs0ln-co_Wx0ln*x_bar0ln)*mh_0ln*d_n); 
    
    
    
    ############b) estimating equation of v #########
    ############ score function of v ##############
    ###### the 2nd part: dl2/dv #######
    sc_v2ln=-v_h0ln/alpha_h0ln;  ### n*1 vector
    
    ##### initial matrix ######
    mat_vln=matrix(0,nrow=D,ncol=n);  
    for (a in 1:n) {
      if(a>1){
        n_0=sum(q[1:(a-1)])+1; ### the start number of cluster a
        n_1=sum(q[1:a]);   ### the end number of cluster a
      }else{
        n_0=1;n_1=q[1];
      }
      mat_vln[,a]=colSums(M_hat0ln[n_0:n_1,]);  #### D*1 vector
    }
    
    ################ the 1st part: dl1/dv ##############
    sc_v1ln=colSums(mat_vln*d_n);  #### n*1 vector
    
    ###############################################################
    ################# final estimate eq of v ################
    dft_vln=sc_v1ln+sc_v2ln; ### n*1 vector
    
    
    ######## estimate A_v in Hessian matrix ############
    ############# the 1st part of -d^2l1/dv^2 ##############
    com_v0ln=v_mat0ln*mh_0ln*d_n; #### D*n matrix
    h_v10ln=-t(v_bar0ln)%*%(com_v0ln);  #### n*n matrix
    ############# the 2nd part of -d^2l1/dv^2 ##############
    hv2_cln=colSums(com_v0ln);  ### n*1 vector
    h_v20ln=diag(hv2_cln);  #### n*n matrix
    ############# the part of -d^2l2/dv^2 ##############
    h_v30ln=diag(1/alpha_h0ln,n,n);  #### n*n matrix
    ############# final estimator of A_v: -d^2l/dv^2 ##############
    Hv_a0ln=h_v10ln+h_v20ln+h_v30ln;  #### n*n matrix
    
    
    ## cross Hessian matrix of beta and v: n*1 vector  ##   
    Mh_bv0ln=(vz_mat0ln-v_mat0ln*x_bar0ln)*mh_0ln*d_n; ### D*n matrix 
    ############# A_bv: -d^2l/dbdv ##############
    Hbv_a0ln=colSums(Mh_bv0ln);  ## n*1 vector 
    
    
    ################### inverse Hessian matrix ###############
    hiv_a0ln=solve(Hv_a0ln); ### the inverse matrix of A_v:n*n matrix
    H1_co0ln=colSums(Hbv_a0ln*hiv_a0ln);  ## n*1 vector
    H2_co0ln=Hb_a0ln-sum(Hbv_a0ln*H1_co0ln);
    #### the 1st part in inverse Hessian matrix ##
    H11_0ln=1/H2_co0ln; 
    #### the 2nd part in inverse Hessian matrix ##
    H12_0ln=-H1_co0ln*H11_0ln; ### n*1 vector
    #### the 4th part in inverse Hessian matrix ##
    H22_0ln=hiv_a0ln+H1_co0ln%*%t(H1_co0ln)/H2_co0ln; ###  log-normal frailty, n*n matrix
    
    
    ######### estimate beta and v together ###########
    beta_h1ln=beta_h0ln+H11_0ln*dft_bln+sum(H12_0ln*dft_vln);
    v_h1ln=v_h0ln+H12_0ln*dft_bln+colSums(H22_0ln*dft_vln);   #### n*1 vector
    
    
    ########### updated basic parts #########
    v_1ln=rep(v_h1ln,q); ### N*1 vector
    cv_1ln=x1*beta_h1ln+v_1ln;    ### eta_hat: N*1 vector
    expeta_h1ln=as.numeric(exp(cv_1ln));  ### exp(eta_hat): N*1 vector
    ############ update common parts (include beta and v) ###############
    Wi_exp1ln=Wi_1km*expeta_h1ln;  #### N*D matrix
    Wx_exp1ln=Wi_x1km*expeta_h1ln;  #### N*D matrix
    co_Wx1ln=colSums(Wx_exp1ln);  #### D*1 vector
    co_Wi1ln=colSums(Wi_exp1ln);  #### D*1 vector
    co_Wxs1ln=colSums(Wi_x1skm*expeta_h1ln);  #### D*1 vector
    
    
    ##### initial matrix: numerator of v_bar #####
    v_mat1ln=matrix(0,nrow=D,ncol=n);  
    vz_mat1ln=matrix(0,nrow=D,ncol=n);  
    for (a in 1:n) {
      if(a>1){
        n_0=sum(q[1:(a-1)])+1; ### the start number of cluster a
        n_1=sum(q[1:a]);   ### the end number of cluster a
      }else{
        n_0=1;n_1=q[1];
      }
      v_mat1ln=matrix(0,nrow=D,ncol=n);  
      vz_mat1ln=matrix(0,nrow=D,ncol=n);  
    }
    
    
    #### update m0_hat, v_bar and x_bar #######  
    if(num_pr0km>0){
      ## since sum(co_Wx0[,wk_0km])=0, we just consider (D-wk_0km) vector
      mh_1ln=(Mc_0km/co_Wi1ln)[-wk_0km]; ## (D-wk_0km)*1 vector 
      v_bar1ln=v_mat1ln[-wk_0km,]/co_Wi1ln[-wk_0km]; ## (D-wk_0km)*n matrix 
      x_bar1ln=(co_Wx1ln/co_Wi1ln)[-wk_0km]; ## (D-wk_0km)*1 vector
      ### final estimator of m0, v_bar and x_bar: let m0[t(D)]=0,v_bar[t(D)]=1/n,x_bar[t(D)]=x_mean
      mh_1ln=c(mh_1ln,mh_wk0); #### D*1 vector 
      v_bar1ln=rbind(v_bar1ln,v_wk0); #### D*n matrix 
      x_bar1ln=c(x_bar1ln,x_wk0); #### D*1 vector
      
    }else{ 
      ### estimator of m0[t(k)], v_bar[t(k)] and x_bar[t(k)]: k=1,...,D ############
      mh_1ln=Mc_0km/co_Wi1ln; ## D*1 vector 
      v_bar1ln=v_mat1ln/co_Wi1ln; ## D*n matrix
      x_bar1ln=co_Wx1ln/co_Wi1ln; ## D*1 vector
    }
    
    
    M_12ln=t(t(Wi_exp1ln)*mh_1ln);
    
    ###### update A_b with frailty in Hessian matrix ######
    Hb_a1ln=sum((co_Wxs1ln-co_Wx1ln*x_bar1ln)*mh_1ln*d_n); 
    
    ## update cross Hessian matrix of beta and v ##  
    Mh_bv1ln=(vz_mat1ln-v_mat1ln*x_bar1ln)*mh_1ln*d_n; ### D*n matrix 
    ############# A_bv: -d^2l/dbdv ##############
    Hbv_a1ln=colSums(Mh_bv1ln);  ## n*1 vector 
    
    
    ######## update first two parts of A_v in Hessian matrix ############
    ############# the 1st part of -d^2l1/dv^2 ##############
    com_v1ln=v_mat1ln*mh_1ln*d_n;
    h_v11ln=-t(v_bar1ln)%*%(com_v1ln);  #### n*n matrix
    ############# the 2nd part of -d^2l1/dv^2 ##############
    hv2_c1ln=colSums(com_v1ln);  ### n*1 vector
    h_v21ln=diag(hv2_c1ln);  #### n*n matrix
    
    
    ########### part of estimate se of beta ############
    
    #### update M(t(k)):N*D matrix #####
    M_hat1ln=M_01km-M_12ln;
    
    
    ##############################################
    ########### estimate Sigma (for beta and v) ############
    
    ####(1) estimate xi_b #######
    ### a) estimate the 1st part of xi_b ####
    xi_b11ln=colSums(t(M_hat1ln*x1)*d_n);  ### N*1 vector
    xi_b12ln=colSums(t(M_hat1ln)*x_bar1ln*d_n);  ### N*1 vector
    ### final estiamtor of 1st part of xi_b ####
    xi_b1ln=xi_b11ln-xi_b12ln;  ### N*1 vector
    ### b) estimate the 2nd part of xi_b ####
    ######## estimator of Q_b ##########
    Qb_kln=colSums(M_cg*xi_b1ln);   ### D_c*1 vector
    ##### final estiamtor of 2nd part of xi_b ######
    xi_b21ln=Qb_kln/n_cg;  ### D_c*1 vector
    xi_b2ln=colSums(xi_b21ln*dMc_hat);  ### N*1 vector
    ######### final estimator of xi_b ###########
    xi_bln=xi_b1ln+xi_b2ln;  ### N*1 vector  
    ######## final estimator of Sigma_b ############
    Hb_sln=sum(xi_bln*xi_bln);  
    
    
    ####(2) estimate xi_v #######
    ####### initial vectors #######
    xi_v1ln=NA;
    for (a in 1:n) {
      if(a>1){
        n_0=sum(q[1:(a-1)])+1; ### the start number of cluster a
        n_1=sum(q[1:a]);   ### the end number of cluster a
      }else{
        n_0=1;n_1=q[1];
      }
      ### a) estimate the 1st part of xi_v ####
      ###### log-normal frailty
      xi_v11ln=colSums(t(M_hat1ln[n_0:n_1,])*d_n);  ### q*1 vector
      xi_v12ln=colSums(t(M_hat1ln[n_0:n_1,])*v_bar1ln[,a]*d_n);  ### q*1 vector
      ### final estiamtor of 1st part of xi_v ####
      xi_v1ln[n_0:n_1]=xi_v11ln-xi_v12ln;  ### q*1 vector
      
    }
    
    ### b) estimate the 2nd part of xi_v ####
    ######## estimator of Q_v: D_c*1 vector ##########
    Qv_k1ln=colSums(M_cg*xi_v1ln);    
    ##### final estiamtor of 2nd part of xi_v ######
    xi_v21ln=Qv_k1ln/n_cg;  
    
    
    ########### estimate alpha by optimizing profile h-likelihood #########
    alpha_h1ln=optimise(halphaln,c(0.00001,1.5))$minimum;  
    
    
    #############################################################
    ############ set convergence conditions ############
    be_h0=c(alpha_h0ln,beta_h0ln,v_h0ln);
    be_h1=c(alpha_h1ln,beta_h1ln,v_h1ln);
    
    num_na=sum(is.na(be_h1)); ### the number of NA in be_h1
    
    if((max(abs(be_h1-be_h0))<=0.00001)|(num_na>0)|(iter>300)){
      div=div+ifelse(num_na>0,1,0)+floor(iter/300);
      break} else{iter=iter+1;
      ###### log-normal frailty
      alpha_h0ln=alpha_h1ln;
      beta_h0ln=beta_h1ln;
      v_h0ln=v_h1ln;
      }
  }    #### iter ends here 
  
  
  ###############################     
  ##################### update inverse Hessian matrix ######################
  ############# the part of -d^2l2/dv^2: n*n matrix ##############
  ### log-normal frailty
  h_v31ln=diag(1/alpha_h1ln,n,n); 
  ############# final estimate of A_v: -d^2l/dv^2 ##############
  Hv_a1ln=h_v11ln+h_v21ln+h_v31ln;  
  
  hiv_a1ln=solve(Hv_a1ln); ### the inverse matrix of A_v
  H1_co1ln=colSums(Hbv_a1ln*hiv_a1ln);  ## n*1 vector
  H2_co1ln=Hb_a1ln-sum(Hbv_a1ln*H1_co1ln);
  #### the 1st part in inverse Hessian matrix ##
  H11_1ln=1/H2_co1ln; 
  
  ############################################################# 
  ## sandwich variance of beta ##
  hiv_b1ln=as.numeric(H11_1ln%*%Hb_sln%*%H11_1ln);
  ##### final estimator of se(beta) #####
  se_betah1ln=sqrt(hiv_b1ln); 
  
  
  ################################################
  ############## outputs ##############
  est_listln=list(beta_hat=beta_h1ln,v_hat=v_h1ln,alpha_hat=alpha_h1ln);
  se_listln=list(se_sw=se_betah1ln);
  
  re_listln=list(estimate=est_listln,se=se_listln);
  
  return(re_listln);
}
