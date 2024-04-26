
##############################################################################
####### R code for MRL models with frailty (penalty method)

### remove the history
rm(list=ls(all=TRUE))

## R packages to use Buckley-James estimate ###
library(Hmisc)
library(lattice)
library(Formula)
library(survival)
library(ggplot2)
library(rms)
## R packages to use hessian function in HQL
library(numDeriv)


### load Newton-Raphson iterative functions ###
source("NRcqlbj.R")
source("NRcqlkm.R")
source("NRhqlln.R")
source("NRnon.R")


#########################################################
############### read data ############
our.data=read.csv('./data.csv',header = TRUE);
attach(our.data);

c_c=names(table(o.id));
nam_c=as.vector(c_c);  ### name of cluster
s_c=as.data.frame(table(o.id))$Freq;
q=as.vector(s_c);  ### size of clusters
n=length(nam_c);  ### number of cluster
N=nrow(our.data);   ### sample size
N_c=ncol(our.data);  ### number of variables in the data
lambda=log(N);   ### tunning parameter

y1=o.surv;    ### observed survival times
del1=e.surv;  ### survival indicator
x1=cov;   ### covariate vector

x_mean=mean(x1);
c_r=1-mean(del1); ## censoring rate


yd=cbind(y1,del1);      ### dataset with survival time and indicator
yd_ord=yd[order(y1),];  ### datasset ordered by observed times 
y_ord=yd_ord[,1]; ### ordered survival time
d_ord=yd_ord[,2]; ### ordered event indicators


######## failure times ##########
t1_ord=y_ord[d_ord==1]; ### ordered failure times
t2_dist=unique(t1_ord); ### distinct failure times
D=length(t2_dist);  ### number of distinct failure times


###  Matrix of observed Set: triangle matrix type
### (1) Matrix of y[ij]>t(k),k=1,2,...,D.
### M_g is the N*D indicator matrix whose (ij,k)th element is 1 if y[ij]>t(k) and 0 otherwise. Here t(k) is the kth smallest distinct failure time, k=1,...,D.
M_g=matrix(0,N,D);
for (i in 1:N){
  M_g[i,]=ifelse(y1[i]>t2_dist,1,0)
}     

### (2) Matrix of  y[ij]==t(k)
### M_e is the N*G indicator matrix whose (ij,k)th element is 1 if y[ij]==t(k) and 0 otherwise.
M_e=matrix(0,N,D);
for (i in 1:N){
  M_e[i,]=ifelse(y1[i]==t2_dist,1,0)
}
d_n=colSums(del1*M_e);  #### D*1 vector: failure numbers at t(k)


##########################################################
######### For censoring time #################
c_ord=y_ord[d_ord==0];  ### ordered censoring times
c_dist=unique(c_ord);  ### the distinct censoring times
D_c=length(c_dist);  #### the number of distinct times 

############# common matrix ########
### (1) Matrix of  y[ij]>=c[l],l=1,2,...,D_c
### M_cg is the N*D_c indicator matrix whose (ij,l)th element is 1 if y[ij]>=c[l] and 0 otherwise.
M_cg=matrix(0,N,D_c);
for (i in 1:N){
  M_cg[i,]=ifelse(y1[i]>=c_dist,1,0)
}
n_cg=colSums(M_cg);  ### censoring times at risk: D_c*1 vector

### (2) Matrix of  y[ij]==c(l),l=1,2,...,D_c.
### M_ce is the N*D_c indicator matrix whose (ij,l)th element is 1 if y[ij]==c(l) and 0 otherwise.
M_ce=matrix(0,N,D_c)
for (i in 1:N){
  M_ce[i,]=ifelse(y1[i]==c_dist,1,0)
}
dN_c=(1-del1)*M_ce;  ### N*D_c matirx
d_cn=colSums(dN_c);  ### censoring number at c(l): D_c*1 vector


###########################################################
################ estimate survival of C ##########
##### (1) Kaplan-Meier method ####
###### estimate censoring rates at c(l): D_c*1 vector
cr_c0=d_cn/n_cg;       
###### D_c*N matrix #######
Mc_hat0=t(M_cg)*cr_c0;  
dMc_hat=t(dN_c)-Mc_hat0; 

### the log-transform 
cr_l=log(1-cr_c0);  ### D_c*1 vector
### sum of cr_l at y[ij]: N*1 vector
cumcr_c=NA; 
for(i in 1:N){
  k_1=which(M_cg[i,]==1)
  cumcr_c[i]=sum(cr_l[k_1])
}

survc_h0km=exp(cumcr_c); ### estimates of G(y[ij]): N*1 vector


##### (2) Buckley-James method #####
## without repeat
### Buckley-James estimate using function 'bj' from 'rms' package ####
dd <- datadist(x1)
options(datadist="dd")
fit=bj(Surv(y1,1-del1)~x1, time.inc=1, control=list(iter.max=200))
res_bj=fit$linear.predictors-log(y1) ##-epsilon

### estimate of sig_eps and coeffecient in AFT model
sigmahat=sd(res_bj[del1==0])*sqrt(6)/pi; ## use the part of data with del=0
coef_aft=coef(fit);


################################
### estimate the distribution parameter of C, derive based on AFT model with BJ estimate ######
###### (1)calculate the errors in AFT model based on BJ estimates
ero_aft=(log(y1)-coef_aft[1]-coef_aft[2]*x1)/sigmahat;  ##for 1 covariate
###### (2) Sort (error,delta) by error, to get error(k).
ero_del=cbind(ero_aft,del1);
ero_del_ord=ero_del[order(ero_aft),];   ### ordered errors
ero_ord=ero_del_ord[,1];    ### observed error 
del_ord=ero_del_ord[,2];    ### corresponding delta
###### (3) Process for finding the distinct error failure times
ero1_ord=ero_ord[del_ord==0];   ###the ordered error failure times
ero1_dist=unique(ero1_ord);  ###the distinct error failure times
D_ero=length(ero1_dist); #### the number of distinct times
###### (4) Matrix of observed Set: triangle matrix type
### (a1) Matrix of ero[i]>ero1[k],k=1,2,...,D_ero
### M_ero_g is the N*D_ero indicator matrix whose (i,k)th element is 1 if ero[i]>_ero1[k] and 0 otherwise. Here ero1[k] is the kth smallest distinct error failure time, k=1,...,D_ero
M_ero_cg=matrix(0,N,D_ero);
for (i in 1:N){
  M_ero_cg[i,]=ifelse(ero_aft[i]>=ero1_dist,1,0)
}     
n_ero=colSums(M_ero_cg);   ##column sums of M1_ero_g: the kth element is the at risk number at ero1[k]

### (a2) Matrix of ero[i]==ero1[k]
### M_e is the N*D_ero indicator matrix whose (i,k)th element is 1 if y1[i]==ero1[k] and 0 otherwise.
M_ero_e=matrix(0,N,D_ero);
for (i in 1:N){
  M_ero_e[i,]=ifelse(ero_aft[i]==ero1_dist,1,0)
}
dN_ero=(1-del1)*M_ero_e; ### N*D_ero matirx
d_ero=colSums(dN_ero);  #### D_ero*1 vector: failure numbers at ero1[k]

### a) estimate censoring rates at ero_dist[l]: D_ero*1 vector
cr_ero=d_ero/n_ero; 
cr_ero_l=log(1-cr_ero);  ### the log-transform 
### sum of cr_ero_l at y1[i]: N*1 vector
cumcr_ero=NA; 
for(i in 1:N){
  k_ero_1=which(M_ero_cg[i,]==1)
  cumcr_ero[i]=sum(cr_ero_l[k_ero_1])
}
### b) estimate censoring survival 
survc_h0=exp(cumcr_ero); ### estimates of G_c(ero[i]): n*1 vector

M_ero_hat0=t(M_ero_cg)*cr_ero;  ### D_ero*N matrix
dM_ero_hat=t(dN_ero)-M_ero_hat0;  ### D_ero*N matrix


######### the inverse weight (IPCW) #############
## (1) KM-based
r_dkm=ifelse(survc_h0km==0,0,del1/survc_h0km);  ### N*1 vector
########## final IPCW ##############
Wi_1km=r_dkm*M_g;  ### N*D matrix
### common conditions ###
in_par0km=colSums(Wi_1km);  #### D*1 vector
wk_0km=which(in_par0km==0);
num_pr0km=length(wk_0km);
x_wk0km=matrix(rep(x_mean,num_pr0km),nrow=num_pr0km,byrow=T); ## num_pr0*p matrix

############ common parts ###############
Wi_ykm=Wi_1km*y1;    #### N*D matrix
Wi_tkm=t(t(Wi_1km)*t2_dist); ### N*D matrix
M_01km=Wi_ykm-Wi_tkm;    #### N*D matrix
Mc_0km=colSums(M_01km);  #### D*1 vector

#### covariates: N*D matrix ######
Wi_x1km=Wi_1km*x1;     
Wi_x1skm=Wi_x1km*x1;    

## (2) BJ-based
r_d=ifelse(survc_h0==0,0,del1/survc_h0);  ### N*1 vector
########## final IPCW ##############
Wi_1=r_d*M_g;  ### N*D matrix

### common conditions ###
in_par0=colSums(Wi_1);  #### D*1 vector
wk_0=which(in_par0==0);
num_pr0=length(wk_0);
x_wk0=matrix(rep(x_mean,num_pr0),nrow=num_pr0,byrow=T); ## num_pr0*p matrix

############ common parts ###############
Wi_y=Wi_1*y1;    #### N*D matrix
Wi_t=t(t(Wi_1)*t2_dist); ### N*D matrix
M_01=Wi_y-Wi_t;    #### N*D matrix
Mc_0=colSums(M_01);  #### D*1 vector

#### covariates: N*D matrix ######
Wi_x1=Wi_1*x1;    
Wi_x1s=Wi_x1*x1;    


##################################################################
################### estimate starts here ###################

##### estimate beta & v and their SEs ####
### initial values
alpha_h0=0.1;
beta_h0=0;
v_h0=rep(0,n);

p=length(beta_h0);

#### without penalty with KM-IPCW #####
re_non=NR_non(beta_h0);
########(i) estimate beta ########
beta_h1non=re_non$beta_hatnon;
########(ii) estimate se(beta) ########
se_betah1non=re_non$se_swnon; ## sandwich

#### penalty ####
#### (1) based on KM estimate ####
re_km=CQL_km(beta_h0,v_h0);
########(i) estimate beta and v ########
beta_h1km=re_km$estimate$beta_hat;
v_h1km=re_km$estimate$v_hat;
########(ii) estimate se(beta) ########
se_betah1km=re_km$se$se_sw;  ### the main part (1st part) of sandwich


#### (2) based on BJ estimate ####
re_bj=CQL_bj(beta_h0,v_h0);
########(i) estimate beta and v ########
beta_h1bj=re_bj$estimate$beta_hat;
v_h1bj=re_bj$estimate$v_hat;
########(ii) estimate se(beta) ########
se_betah1bj=re_bj$se$se_sw;  ### the main part (1st part) of sandwich


#### (3) based on HQL estimate with lognormal ####
re_hql=HQL_lnkm(beta_h0,v_h0,alpha_h0);  
########(i) estimate beta and v ########
beta_h1hqlln=re_hql$estimate$beta_hat;
v_h1hqlln=re_hql$estimate$v_hat;
########(ii) estimate se(beta) ########
se_betah1hqlln=re_hql$se$se_sw;  ### the main part (1st part) of sandwich


#### store the results #######
cr0=round(c_r,4); ## censoring rate

## without penalty
beta_hatnon=round(beta_h1non,4);
beta_se1non=round(se_betah1non,4);
renon=cbind(n,N,cr0,beta_hatnon,beta_se1non);

## with penalty
#### (1) based on KM estimate ####
beta_hatkm=round(beta_h1km,4);
beta_se1km=round(se_betah1km,4);
rekm=cbind(n,N,cr0,beta_hatkm,beta_se1km);

#### (2) based on BJ estimate ####
beta_hatbj=round(beta_h1bj,4);
beta_se1bj=round(se_betah1bj,4);
rebj=cbind(n,N,cr0,beta_hatbj,beta_se1bj);

#### (3) based on HQL estimate with lognormal ####
beta_hathqlln=round(beta_h1hqlln,4);
beta_se1hqlln=round(se_betah1hqlln,4);
rehqlln=cbind(n,N,cr0,beta_hathqlln,beta_se1hqlln);


#########################################################
### return results
rebj;rekm;rehqlln;renon;

##################### R code is over ####################
