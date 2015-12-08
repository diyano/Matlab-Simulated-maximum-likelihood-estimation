function [prob_a, stage2_prob_early_ifa,stage2_prob_early_ifb, eu_a, eu_b, eu_a_withflex, eu_b_withflex] ...
    = solvemodel(theta_beta,theta_delta,theta_sigma,theta_omega,phi,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws,persnr_cont,sophisticated)

% Reorganize the data. 

N=size(unique(persnr_cont),1);
Nobs=size(a1,1) ;

Ndraws=size(X_beta,3);

theta_beta_array = repmat(reshape(theta_beta,1,[],1),Nobs,1,Ndraws);
theta_delta_array = repmat(reshape(theta_delta,1,[],1),Nobs,1,Ndraws);
theta_sigma_array = repmat(reshape(theta_sigma,1,[],1),Nobs,1,Ndraws);
theta_omega_array = repmat(reshape(theta_omega,1,[],1),Nobs,1,Ndraws);

R=[phi(1) 0 0 0;phi(2) phi(3) 0 0; phi(4) phi(5) phi(6) 0; phi(7) phi(8) phi(9) phi(10)];
% vcov = R'*R
numberrc=size(R,1);
% compare the following to vcov
% cov(standardnormaldraws*R)

normaldraws_temp=reshape((standardnormaldraws*R)',numberrc,N,[]);
normaldraws_permuted=permute(normaldraws_temp,[2 1 3]);
clear temp
% cov(squeeze(draws(7,:,:))')

normaldraws_obs=normaldraws_permuted(persnr_cont,:,:);

beta_matrix=squeeze(sum(X_beta.*theta_beta_array,2) + normaldraws_obs(:,1,:));
delta_matrix=max(0,1-squeeze(sum(X_delta.*theta_delta_array,2) + normaldraws_obs(:,2,:)).^2);
sigma_e_matrix=max(1e-1,min(squeeze((sum(X_sigma.*theta_sigma_array,2) + normaldraws_obs(:,3,:)).^2),10));
omega_matrix=squeeze(sum(X_omega.*theta_omega_array,2) + normaldraws_obs(:,4,:));


%% choice in stage 2

% from the perspective of stage 2, so when in period 2, hence beta*delta

% utility early minus late in period 1 if a chosen
diff_a = repmat(a1,1,Ndraws)-beta_matrix.*delta_matrix.*repmat(a2,1,Ndraws);

% utility early minus late in period 1 if b chosen
diff_b = repmat(b1,1,Ndraws)-beta_matrix.*delta_matrix.*repmat(b2,1,Ndraws);
diff_b(typeq==4,:)=-Inf;
diff_b(typeq==5,:)=Inf;

stage2_prob_early_ifa = cdf('logistic',diff_a./sigma_e_matrix,0,1);
stage2_prob_early_ifb = cdf('logistic',diff_b./sigma_e_matrix,0,1);


%% expected utilities
if sophisticated==1
eu_a_logsum = beta_matrix.*delta_matrix.*...
    (sigma_e_matrix .*exp(1) + sigma_e_matrix .* log(exp(repmat(a1,1,Ndraws)./sigma_e_matrix)...
    +exp(beta_matrix.*delta_matrix.*repmat(a2,1,Ndraws)./sigma_e_matrix))) + ...
    beta_matrix.*delta_matrix.*delta_matrix.*(1-beta_matrix).*...
    exp(beta_matrix.*delta_matrix.*repmat(a2,1,Ndraws)./sigma_e_matrix) ./ (exp(repmat(a1,1,Ndraws)./sigma_e_matrix)+exp(beta_matrix.*delta_matrix.*repmat(a2,1,Ndraws)./sigma_e_matrix));

eu_b_logsum = beta_matrix.*delta_matrix.*...
    (sigma_e_matrix .*exp(1) + sigma_e_matrix .* log(exp(repmat(b1,1,Ndraws)./sigma_e_matrix)...
    +exp(beta_matrix.*delta_matrix.*repmat(b2,1,Ndraws)./sigma_e_matrix)))+ ...
    beta_matrix.*delta_matrix.*delta_matrix.*(1-beta_matrix).*...
    exp(beta_matrix.*delta_matrix.*repmat(b2,1,Ndraws)./sigma_e_matrix) ./ (exp(repmat(b1,1,Ndraws)./sigma_e_matrix)+exp(beta_matrix.*delta_matrix.*repmat(b2,1,Ndraws)./sigma_e_matrix));
end
if sophisticated==0
eu_a_logsum = beta_matrix.*delta_matrix.*...
    (sigma_e_matrix .*exp(1) + sigma_e_matrix .* log(exp(repmat(a1,1,Ndraws)./sigma_e_matrix)...
    +exp(delta_matrix.*repmat(a2,1,Ndraws)./sigma_e_matrix)));

eu_b_logsum = beta_matrix.*delta_matrix.*...
    (sigma_e_matrix .*exp(1) + sigma_e_matrix .* log(exp(repmat(b1,1,Ndraws)./sigma_e_matrix)...
    +exp(delta_matrix.*repmat(b2,1,Ndraws)./sigma_e_matrix)));
end
typeq_array = repmat(typeq,1,Ndraws);

payoff_a_early_type1=repmat(a1,1,Ndraws);
payoff_b_late_type1=beta_matrix.*delta_matrix.*repmat(b1,1,Ndraws);

payoff_a_early_type2=beta_matrix.*delta_matrix.*repmat(a1,1,Ndraws);
payoff_b_late_type2=beta_matrix.*delta_matrix.*delta_matrix.*repmat(b1,1,Ndraws);

payoff_b_late_type4=beta_matrix.*delta_matrix.*delta_matrix.*repmat(b2,1,Ndraws);
payoff_b_early_type5=beta_matrix.*delta_matrix.*repmat(b1,1,Ndraws);

eu_a=NaN(Nobs,Ndraws);
eu_a(typeq_array==1)=payoff_a_early_type1(typeq_array==1);
eu_a(typeq_array==2)=payoff_a_early_type2(typeq_array==2);
eu_a(typeq_array==4|typeq_array==5|typeq_array==6)=...
eu_a_logsum(typeq_array==4|typeq_array==5|typeq_array==6);

eu_b=NaN(Nobs,Ndraws);
eu_b(typeq_array==1)=payoff_b_late_type1(typeq_array==1);
eu_b(typeq_array==2)=payoff_b_late_type2(typeq_array==2);
eu_b(typeq_array==4)=payoff_b_late_type4(typeq_array==4);
eu_b(typeq_array==5)=payoff_b_early_type5(typeq_array==5);
eu_b(typeq_array==6)=eu_b_logsum(typeq_array==6);
%% difference
flexible_a = repmat(typeq==4|typeq==5|typeq==6,1,Ndraws);
flexible_b = repmat(typeq==6,1,Ndraws);
%% EU with flex
eu_a_withflex= eu_a + flexible_a .* omega_matrix;
eu_b_withflex= eu_b + flexible_b .* omega_matrix;
%% probability a

prob_a = cdf('logistic',(eu_a+flexible_a.*omega_matrix-eu_b-flexible_b.*omega_matrix)./sigma_e_matrix);

%{
left=eu_a+flexible_a.*omega_matrix;
right=eu_b-flexible_b.*omega_matrix;
temp=[prob_a(:,1) left(:,1) right(:,1) sigma_e_matrix(:,1) eu_a(:,1) flexible_a(:,1) omega_matrix(:,1) eu_b(:,1) flexible_b(:,1) omega_matrix(:,1) ];
%}
lastline=0;