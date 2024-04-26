
#####################################################
##### Newton-Raphson iteration: model without penalty ######
#### estimate beta and se(beta) based on KM estiamtor ###

NR_non=function(betanon){
  
  ########################################################
  ##### PART 1: estimate beta #######
  ######### iteration starts here  ##########
  iter=1;div=0;   ### iteration and diverging times
  repeat
  { 
    #############################
    #####(i) basic parts ####
    c_0non=x1*betanon;   ### N*1 vector
    expeta_h0non=as.numeric(exp(c_0non));  ### exp(X*beta): N*1 vector
    ############a) common parts including beta ###############
    Wi_exp0non=Wi_1km*expeta_h0non;  #### N*D matrix
    Wx_exp0non=Wi_x1km*expeta_h0non;  #### N*D matrix
    co_Wi0non=colSums(Wi_exp0non);  #### D*1 vector
    co_Wx0non=colSums(Wx_exp0non);  #### D*1 vector
    co_Wxs0non=colSums(Wi_x1skm*expeta_h0non);  #### D*1 vector
    ######b) estimate m0[t(k)] and x_bar[t(k)]: k=1,...,D ############
    mh_0non=Mc_0km/co_Wi0non;
    x_bar0non=co_Wx0non/co_Wi0non;
    if(num_pr0km>0){
      mh_0non[wk_0km]=0;
      x_bar0non[wk_0km]=x_mean;
    }
    ####c) estimate M(t(k)): N*D matrix #####
    M_02non=t(t(Wi_exp0non)*mh_0non);
    M_hat0non=M_01km-M_02non;
    
    
    
    #####################################
    #####(ii) estimate beta #######
    
    ####### common matrix ##########
    mhd_0=mh_0non*d_n; #### D*1 vector
    
    ############a) estimating equation of beta #########
    sc_bnon=colSums(M_hat0non*x1); #### D*1 vector
    dft_bnon=sum(sc_bnon*d_n); 
    ####b) estimate Hessian matrix: Hbb ########
    Hb_a0non=sum((co_Wxs0non-co_Wx0non*x_bar0non)*mhd_0); 
    #######c) final estimate of beta ###########
    beta_h1non=betanon+dft_bnon/Hb_a0non;
    
    
    
    ###########################################
    ############(iii) set convergence conditions ############
    num_na=sum(is.na(beta_h1non)); ### the number of NA in theta_h1
    
    if((abs(beta_h1non-betanon)<=0.00001)|(num_na>0)|(iter>300)){
      div=div+ifelse(num_na>0,1,0)+floor(iter/300);
      break} else{iter=iter+1;
      betanon=beta_h1non;
      }
  }    #### iter ends here 
  
  
  
  ##################################################################     
  ########### PART 2: estimate se of beta ############
  
  ###############################################
  ###########(ii) estimate Sigma_b in covariance ############
  
  #####a) estimate xi_b ####
  #####a1) estimate the 1st part of xi_b ####
  xi_b11non=colSums(t(M_hat0non*x1)*d_n);  ### N*1 vector
  xi_b12non=colSums(t(M_hat0non)*x_bar0non*d_n);  ### N*1 vector
  ### final estiamtor of 1st part of xi_b ####
  xi_b1non=xi_b11non-xi_b12non;  ### N*1 vector
  ######a2) estimate the 2nd part of xi_b ####
  ######## estimates of Q_b ##########
  Qb_knon=colSums(M_cg*xi_b1non);   ### D_c*1 vector
  ##### final estiamtor of 2nd part of xi_b ######
  xi_b21non=Qb_knon/n_cg;  ### D_c*1 vector
  xi_b2non=colSums(xi_b21non*dMc_hat);  ### N*1 vector
  ######### final estimates of xi_b ###########
  xi_bnon=(xi_b1non+xi_b2non);  ### N*1 vector 
  
  ######## final estimate of Sigma_b ############
  Hb_snon=sum(xi_bnon*xi_bnon);
  
  
  ###########################################
  ######(iii) update inverse Hessian matrix ######
  ###a) update Hbb in Hessian matrix ######
  Hb_a10non=co_Wxs0non-(co_Wx0non*x_bar0non);
  Hb_a1non=sum(Hb_a10non*mhd_0);
  #######b) inverse Hessian matrix ######
  H11_1non=1/Hb_a1non; #### Hbb^inv
  
  
  ###############################################
  ######(iv) estimate se(beta) ####
  se_betah1non=H11_1non*sqrt(Hb_snon); ## sandwich
  
  
  ############################
  #### outputs ####
  rnon_list=list(beta_hatnon=beta_h1non,se_swnon=se_betah1non);
  return(rnon_list);
}
