function [nll,ngrad,nHess,nobs,N_stage1_prob_a,S_stage1_prob_a,stage2_prob_early_ifa,stage2_prob_early_ifb,N_eu_a,N_eu_b,S_eu_a,S_eu_b] ...
    = negloglik(theta_beta1,theta_delta1,theta_sigma1,theta_omega1,phi1,theta_beta2,theta_delta2,theta_sigma2,theta_omega2,phi2,...
    pi,choose_a,choose_a_stage2,a1,a2,b1,b2,typeq,X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1, standardnormaldraws2, persnr_cont)

Ndraws = size(X_delta,3);

% with or without time preference elicitation
%whichobs = typeq>=4;
whichobs = typeq>0;

% calculate number observations per person
[~,persnr_firstobs,~] = unique(persnr_cont(whichobs));
obs_per_person=persnr_firstobs(2)-persnr_firstobs(1);

% solve model for naive
theta_beta=theta_beta1;
theta_delta=theta_delta1;
theta_sigma=theta_sigma1;
theta_omega=theta_omega1;
phi=phi1;
standardnormaldraws=standardnormaldraws1;
sophisticated=0;

[N_stage1_prob_a,stage2_prob_early_ifa,stage2_prob_early_ifb,N_eu_a,N_eu_b] ...
    = solvemodel(theta_beta,theta_delta,theta_sigma,theta_omega,phi,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws,persnr_cont,sophisticated);

% solve model for sophisticated
theta_beta=theta_beta2;
theta_delta=theta_delta2;
theta_sigma=theta_sigma2;
theta_omega=theta_omega2;
phi=phi2;
standardnormaldraws=standardnormaldraws2;
sophisticated=1;

[S_stage1_prob_a,stage2_prob_early_ifa,stage2_prob_early_ifb,S_eu_a,S_eu_b] ...
    = solvemodel(theta_beta,theta_delta,theta_sigma,theta_omega,phi,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws,persnr_cont,sophisticated);

    % choose_a : whether you choose a or not in the first stage (2709 x 1)  
    % choose_a_expanded: whether you choose a or not in the first stage (2709x1x100)

choose_a_expanded = repmat(choose_a,1,Ndraws);    
    
% reshape choice data first stage
    % prob_choice: first stage likelihood contribution 2709x1x100
    % prob_choice_expanded: first stage likelihood contribution 21 x 129 x 100
N_prob_choice = choose_a_expanded.*N_stage1_prob_a+(1-choose_a_expanded).*(1-N_stage1_prob_a);
N_prob_choice_reshaped = reshape(N_prob_choice(whichobs,:),obs_per_person,[],Ndraws);

S_prob_choice = choose_a_expanded.*S_stage1_prob_a+(1-choose_a_expanded).*(1-S_stage1_prob_a);
S_prob_choice_reshaped = reshape(S_prob_choice(whichobs,:),obs_per_person,[],Ndraws);
% S_prob_choice_reshaped has zeros
% temp=S_prob_choice_reshaped(:,:,1)
% S_prob_choice has zeros

%{
onemin_choose_a_expanded=1-choose_a_expanded;
onemin_S_stage1_prob_a=1-S_stage1_prob_a;
temp=[typeq S_prob_choice(:,1) choose_a_expanded(:,1) S_stage1_prob_a(:,1) onemin_choose_a_expanded(:,1) onemin_S_stage1_prob_a(:,1)]
%row 1250 has zero for example. 
%S_stage1_prob_a cannot be exactly 1 or 0
%}

% reshape choice data second stage
    % the answer given in second stage (1 late 0 early n/a no answer):
    % choose_a_stage2: 2709 x1
    % choose_a_expanded_stage2: 2709 x 1 x100
choose_a_expanded_stage2 = repmat(choose_a_stage2,1,Ndraws);

    % prob_a_stage2 = prob early in second stage (taking chosen first stage option into account)

prob_a_stage2 = choose_a_expanded.*stage2_prob_early_ifa+(1-choose_a_expanded).*stage2_prob_early_ifb;
    
    % prob_choice_stage2: likelihood contribution from the second stage choice: 2709 x 1 x 100
prob_choice_stage2 = choose_a_expanded_stage2.*prob_a_stage2+(1-choose_a_expanded_stage2).*(1-prob_a_stage2);
    
    % prob_choice_reshaped_stage2: likelihood contribution from the second stage choice = 21 x 129 x 100

 prob_choice_reshaped_stage2 = reshape(prob_choice_stage2(whichobs,:),obs_per_person,[],Ndraws);
    
    % temp: likelihood contribution from the second stage choice (n/a replaced by 0)= 21 x 129 x 100

temp = prob_choice_reshaped_stage2;
temp(isnan(temp)) = 0;

    % prob_choice_reshaped_stage2_pick: summing over observations gives
    % likelihood contribution for the second stage per individual: 1 x 129 x 100 
prob_choice_reshaped_stage2_pick = squeeze(sum(temp,1));

% product over time periods
    
    % product over time periods for stage 1
    % prob_choice_reshaped: first stage likelihood contribution per question 21 x 129 x 100
    % prob_allchoice_stage1: first stage likelihood contribution, all questions 129 x 100
N_prob_allchoice_stage1 = squeeze(prod(N_prob_choice_reshaped,1));
S_prob_allchoice_stage1 = squeeze(prod(S_prob_choice_reshaped,1));
% there are zeros in N_prob_allchoice_stage1

    % nobs_stage1. number of subjects x number of draws= 129x100 
nobs_stage1=size(N_prob_choice_reshaped,1)*size(N_prob_choice_reshaped,2);

% stage 2
    % prob_choice_reshaped_stage2: likelihood contribution from the second stage choice = 21 x 129 x 100
    % makes_choice_stage2: if, in the second stage, likelihood contribution
    % is not all empty, e makes a choice in the second stage (makes_choice_stage2=1)
    % isnan term is 21x129x1. summation is taken over rows leading to 1x129x1
    % after transpose makes_choice_stage2= 129x1x1
    % makes_choice_stage2_expanded: 129x1x100 (1 if makes a second stage choice)
makes_choice_stage2 = sum(1-isnan(prob_choice_reshaped_stage2(:,:,1)))';
makes_choice_stage2_expanded = repmat(makes_choice_stage2,1,Ndraws);

    % nobs_stage2 : summation over 129 subjects gives the number of responses in the second stage
    % 1 x 1 x 1
nobs_stage2=sum(makes_choice_stage2);

    % prob_allchoice_stage1: first stage likelihood contribution, per subject 129 x 100
    % makes_choice_stage2_expanded: 129x1x100 (1 if makes a second stage choice)
    % prob_choice_reshaped_stage2_pick: 
    % second stage likelihood contribution per subject: 129 x 100
    
N_prob_allchoice = N_prob_allchoice_stage1.*( (1-makes_choice_stage2_expanded) ...
    + makes_choice_stage2_expanded.*prob_choice_reshaped_stage2_pick);

S_prob_allchoice = S_prob_allchoice_stage1.*( (1-makes_choice_stage2_expanded) ...
    + makes_choice_stage2_expanded.*prob_choice_reshaped_stage2_pick);

weighted_prob_allchoice= S_prob_allchoice .* pi + N_prob_allchoice .* (1-pi);
% average over random draws, then take log, then sum over individuals
nll = - sum(log(mean(weighted_prob_allchoice,2)));

% additional output arguments
ngrad=[];
nHess=[];
nobs=nobs_stage1+nobs_stage2;



