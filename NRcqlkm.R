
#####################################################
##### Newton-Raphson iteration: conditional inference based on KM ######

CQL_km=function(beta_h0,v_h0){
  
  ########################################################
  ##### PART 1: estimate beta and v #######
  
  ######### iteration starts here  ##########
  iter=1;div=0;   ### iteration and diverging times
  repeat
  { 
    #############################
    #####(i) basic parts ####
    v_0=rep(v_h0,each=q); ### N*1 vector
    cv_0=x1*beta_h0+v_0;    ### eta_hat: N*1 vector
    expeta_h0=as.numeric(exp(cv_0));  ### exp(eta_hat): N*1 vector
    ###########a) common parts include beta and v #########
    Wi_exp0=Wi_1km*expeta_h0;  #### N*D matrix
    Wx_exp0=Wi_x1km*expeta_h0;  #### N*D matrix
    co_Wi0=colSums(Wi_exp0);  #### D*1 vector
    co_Wx0=colSums(Wx_exp0);  #### D*1 vector
    co_Wxs0=colSums(Wi_x1skm*expeta_h0);  #### D*1 vector
    #####b) initial matrix: numerator of v_bar #####
    v_mat0=matrix(0,nrow=D,ncol=n);  
    vx_mat0=matrix(0,nrow=D,ncol=n);  
    for (a in 1:n) {
      n_0=(a-1)*q+1;     ### the start number of cluster a
      n_1=a*q;           ### the end number of cluster a
      v_mat0[,a]=colSums(Wi_exp0[n_0:n_1,]);  #### D*1 vector
      vx_mat0[,a]=colSums(Wx_exp0[n_0:n_1,]);  #### D*1 vector  
    }
    
    ######c) estimate m0[t(k)], x_bar[t(k)] and v_bar[t(k),n]: k=1,...,D ############
    mh_0=Mc_0km/co_Wi0;
    x_bar0=co_Wx0/co_Wi0;
    v_bar0=v_mat0/co_Wi0; ## D*n matrix 
    if(num_pr0km>0){
      mh_0[wk_0km]=0;
      x_bar0[wk_0km]=x_mean;
      v_bar0[wk_0km,]=1/n;
    }
    ####d) estimate M(t(k)): N*D matrix #####
    M_02=t(t(Wi_exp0)*mh_0);
    M_hat0=M_01km-M_02;
    sc_com=colSums(t(M_hat0)*d_n); #### N*1 vector
    
    
    
    #####################################
    #####(ii) estimating equations of beta and v #######
    
    ############a) estimating eq of beta #########
    dft_b=sum(sc_com*x1); 
    ############b) estimating eq of v #########
    ############ penalty is sum(v^2/2) ##############
    sc_v2=v_h0;  ### n*1 vector
    ################ the 1st part: dl1/dv ##############
    sc_v1=NA;   ##### initial vector
    for (a in 1:n) {
      n_0=(a-1)*q+1;     ### the start number of cluster a
      n_1=a*q;           ### the end number of cluster a
      sc_v1[a]=sum(sc_com[n_0:n_1]); 
    }
    ################# final estimate eq of v ################
    dft_v=sc_v1-lambda*sc_v2; ### n*1 vector
    
    
    ############################################################
    ####(iii) estimate inverse Hessian matrixes ########
    
    ####### common matrix ##########
    mhd_0=mh_0*d_n; #### D*1 vector
    
    #######a) estimate Hessian matrix: Hbb ############
    Hb_a0=sum((co_Wxs0-co_Wx0*x_bar0)*mhd_0); 
    
    ########b) estimate Hessian matrix: Hvv ############
    ############# the 1st part of -d^2l1/dv^2 ##############
    com_v0=v_mat0*mhd_0; #### D*n matrix
    hv1_c=colSums(com_v0);  ### n*1 vector
    h_v10=diag(hv1_c);  #### n*n matrix
    ############# the 2nd part of -d^2l1/dv^2 ##############
    h_v20=t(v_bar0)%*%(com_v0);  #### n*n matrix
    ############ penalty is sum(v^2/2) ##############
    h_v30=diag(1,n,n);  #### n*n matrix
    
    ############# final estimates of Hvv: -d^2l/dv^2 ##############
    Hv_a0=h_v10-h_v20+lambda*h_v30;  #### n*n matrix
    
    #######c) cross Hessian matrix: Hbv ########   
    Mh_bv0=(vx_mat0-v_mat0*x_bar0)*mhd_0; ### D*n matrix 
    ############# Hbv: -d^2l/dbdv ##############
    Hbv_a0=colSums(Mh_bv0);  ## n*1 vector 
    
    #######d) inverse Hessian matrix #########
    hiv_av0=solve(Hv_a0); ### the inverse matrix of Hvv:n*n matrix
    Hbv_co0=colSums(Hbv_a0*hiv_av0);  ## n*1 vector
    Hbb_co0=Hb_a0-sum(Hbv_a0*Hbv_co0);
    Hvv_co0=Hv_a0-Hbv_a0%*%t(Hbv_a0)/Hb_a0;
    ####d1) Hbb_inv: the 1st part in inverse Hessian matrix ##
    H11_0=1/Hbb_co0; 
    ####d2) Hbv_inv: the 2nd part in inverse Hessian matrix ##
    H21_0=-Hbv_co0*H11_0; ### n*1 vector
    ####d3) Hvv_inv: the 4th part in inverse Hessian matrix ##
    H22_0=solve(Hvv_co0);
    
    
    ###############################################
    #############(iv) final estimates ##########
    #########a) estimate beta and v together ###########
    beta_h1=beta_h0+H11_0*dft_b+sum(H21_0*dft_v);
    v_h1=v_h0+H21_0*dft_b+colSums(H22_0*dft_v);   #### n*1 vector
    
    ############b) set convergence conditions ############
    theta_h0=c(beta_h0,v_h0);
    theta_h1=c(beta_h1,v_h1);
    
    num_na=sum(is.na(theta_h1)); ### the number of NA in theta_h1
    
    if((max(abs(theta_h1-theta_h0))<=0.00001)|(num_na>0)|(iter>300)){
      div=div+ifelse(num_na>0,1,0)+floor(iter/300);
      break} else{iter=iter+1;
      beta_h0=beta_h1;
      v_h0=v_h1;
      }
  }    #### iter ends here 
  
  
  ##################################################################     
  ########### PART 2: estimate se of beta ############
  
  ###############################################
  ###########(i) estimate Sigma in covariance ############
  #####a) estimate Sigma_b
  ####a1) estimate the 1st part of xi_b ####
  xi_b11=colSums(t(M_hat0*x1)*d_n);  ### N*1 vector
  xi_b12=colSums((t(M_hat0)*x_bar0)*d_n);  ### N*1 vector
  ### final estiamtor of 1st part of xi_b ####
  xi_b1=xi_b11-xi_b12;  ### N*1 vector
  ####a2) estimate the 2nd part of xi_b ####
  ######## estimate of Q_b ##########
  Qb_k=colSums(M_cg*xi_b1);   ### D_c*1 vector
  ##### final estiamtor of 2nd part of xi_b ######
  xi_b21=Qb_k/n_cg;  ### D_c*1 vector
  xi_b2=colSums(xi_b21*dMc_hat);  ### N*1 vector
  ######### final estimates of xi_b ###########
  xi_b=(xi_b1+xi_b2);  ### N*1 vector
  
  ######## final estimates of Sigma_b ############
  Hb_s=sum(xi_b*xi_b);  
  
  
  ###############################################
  ######(ii) estimate se(beta) #####
  ####a) sandwich variance of beta ######
  hiv_b1=(H11_0*Hb_s)*H11_0;
  ##### final estimates of se(beta) #####
  se_betah1=sqrt(hiv_b1);  ### the main part (1st part) of sandwich
  
  
  ################################################
  ############## outputs ##############
  est_list=list(beta_hat=beta_h1,v_hat=v_h1);
  se_list=list(se_sw=se_betah1);
  r_list=list(num_iter=iter,num_div=div,estimate=est_list,se=se_list);
  return(r_list);
}
