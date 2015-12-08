clear all
close all
load ../input/long.mat

Hess=0; %Hessian on/off
fontsize=15; % font size for figures
[~,firstobs,persnr_cont] = unique(persnr); % balanced panel
persnr = persnr_cont; 
N=size(firstobs,1);
nobs = size(persnr,1);
%% random draws
numberdraws=100;
numberrc=4;
H1 = haltonset(N*numberrc,'Skip',1e3,'Leap',1e2);
H2 = haltonset(N*numberrc,'Skip',1e3,'Leap',1e2);
H1 = scramble(H1,'RR2');
H2 = scramble(H2,'RR2');
Hdraws1=reshape(net(H1,numberdraws),[],numberrc,1);
Hdraws2=reshape(net(H2,numberdraws),[],numberrc,1);
clear H1 H2 numberrc
standardnormaldraws1=norminv(Hdraws1);
standardnormaldraws2=norminv(Hdraws2);
%% specification
X_beta = [ones(nobs,1,numberdraws)];
X_delta = [ones(nobs,1,numberdraws) repmat((beta1+beta2)/2,1,1,numberdraws)];
X_sigma = [ones(nobs,1,numberdraws)];
X_omega = [ones(nobs,1,numberdraws)];

objfun = @(par)negloglik(par(1),[par(2);0],par(3),par(4),[par(5);par(6);par(7);0;0;0;0;0;0;0],...
par(1),[par(2);0],par(3),par(4),[par(5);par(6);par(7);0;0;0;0;0;0;0],par(8),...    
choose_a,choose_a_stage2,a1,a2,b1,b2,typeq,X_beta,X_delta,X_sigma,X_omega,...
    standardnormaldraws1, standardnormaldraws2, persnr_cont);


startvalues=[ 0.14; 0.31;  1.54 ; 0 ;...
   0.09 ; 0.03; 1.2 ;  ...
   0.25]; % 3.208486693451550e+03
%randomstartingvalues
lb = [-Inf(7,1);0];
ub = [Inf(7,1);1];
%% minimization
%initialize
bhat_knitro = [];
nll_knitro = [];

options_fmincon = optimset('Algorithm','interior-point','Display','iter',...
    'GradObj','off','MaxFunEvals',1e4,'AlwaysHonorConstraints','none',...
    'TolX',1e-8,'TolFun',1e-10);

exist nH var;
nH_exists = ans;
if nH_exists==0
[bhat_fmincon,nll_fmincon,~,~,~,~,nH]...
    = fmincon(objfun,startvalues,[],[],[],[],lb,ub,[],options_fmincon);
        if Hess==1
        nH=hessian(objfun, bhat_fmincon);
        end
else
[bhat_fmincon,nll_fmincon,~,~,~,~,~]...
    = fmincon(objfun,startvalues,[],[],[],[],lb,ub,[],options_fmincon);
end

vcov_fmincon = inv(nH/nobs)/nobs;
se_fmincon = sqrt(diag(vcov_fmincon));

%% the figures
% change this part accordingly when you change the model
beep
bhat=bhat_fmincon;
se=se_fmincon;

bhat1=[bhat(1);bhat(2);bhat(3);bhat(4)];
bhat2=[bhat(1);bhat(2);bhat(3);bhat(4)];
se_bhat1=[se(1);se(2);se(3);se(4)];
se_bhat2=[se(1);se(2);se(3);se(4)];

phi1=[bhat(5);bhat(6);bhat(7);0;0;0;0;0;0;0];
phi2=[bhat(5);bhat(6);bhat(7);0;0;0;0;0;0;0];
se_phi1=[se(5);se(6);se(7);0;0;0;0;0;0;0];
se_phi2=[se(5);se(6);se(7);0;0;0;0;0;0;0];

R1=[phi1(1) 0 0 0;phi1(2) phi1(3) 0 0; phi1(4) phi1(5) phi1(6) 0; phi1(7) phi1(8) phi1(9) phi1(10)];
R2=[phi2(1) 0 0 0;phi2(2) phi2(3) 0 0; phi2(4) phi2(5) phi2(6) 0; phi2(7) phi2(8) phi2(9) phi2(10)];

rc_vcov1 = R1'*R1;
rc_vcov2 = R2'*R2;

randn('seed',1039)
par_draws1=mvnrnd(bhat1',rc_vcov1,1000);
randn('seed',1039)
par_draws2=mvnrnd(bhat2',rc_vcov2,1000);
%% figures and tables
beep
% multistart
fig2 % distribution of random coefficients
fig3 % 21 questions future choice task
fig4 % second stage predictions
fig8 % beta delta relationship
tab1 % estimation results
tab1b % estimation results transformed
tab2 % removing omega and sigma

%% report the results
beep
[bhat_fmincon se_fmincon]
[nll_fmincon]
