% % Categories: Naive(N)/Sophisticated(S) , With Omega (O)/ Without Omega
% (o) , With Sigma (S)/ Without Sigma (s)
close all

%% setting sigma and omega to 0 where necessary
    phi1OS=phi1;
    phi2OS=phi2;
    phi1Os=phi1;
    phi2Os=phi2;
    phi1oS=phi1;
    phi2oS=phi2;
    phi1os=phi1;
    phi2os=phi2;

    phi1Os(4:6)=0;
    phi2Os(4:6)=0;
    phi1os(4:6)=0;
    phi2os(4:6)=0;
    phi1oS(7:10)=0;
    phi2oS(7:10)=0;
    phi1os(7:10)=0;
    phi2os(7:10)=0;
    
    bhat1OS=bhat1;
    bhat2OS=bhat2;
    bhat1Os=bhat1;
    bhat2Os=bhat2;
    bhat1oS=bhat1;
    bhat2oS=bhat2;
    bhat1os=bhat1;
    bhat2os=bhat2;

    bhat1Os(3)=0;
    bhat2Os(3)=0;
    bhat1os(3)=0;
    bhat2os(3)=0;
    bhat1oS(4)=0;
    bhat2oS(4)=0;
    bhat1os(4)=0;
    bhat2os(4)=0;
%
draws_beta1 = max(0,1-par_draws1(:,1).^2);
draws_delta1 = max(0,1-par_draws1(:,2).^2);
draws_sigma1 = max(1e-1,min(par_draws1(:,3).^2,10));
draws_omega1 = par_draws1(:,4);

%example. this can be commented out
a1=10.*ones(size(a1,1),size(a1,2));
a2=a1./mean(draws_delta1,1);

a1_=-Inf(size(a1,1),size(a1,2));
a2_=-Inf(size(a2,1),size(a2,2));

%% Naive
sophisticated=0;

    % Naive(N), With Omega(O), With Sigma(S)
    [~,~,~,~,~,NOSA_eu_a,NOSA_eu_b] ...
    = solvemodel(bhat1OS(1),[bhat1OS(2);0],bhat1OS(3),bhat1OS(4),phi1OS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NOSB_eu_a,NOSB_eu_b] ...
    = solvemodel(bhat1OS(1),[bhat1OS(2);0],bhat1OS(3),bhat1OS(4),phi1OS,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NOSC_eu_a,NOSC_eu_b] ...
    = solvemodel(bhat1OS(1),[bhat1OS(2);0],bhat1OS(3),bhat1OS(4),phi1OS,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), With Omega(O), Without Sigma(s)
    [~,~,~,~,~,NOsA_eu_a,NOsA_eu_b] ...
    = solvemodel(bhat1Os(1),[bhat1Os(2);0],bhat1Os(3),bhat1Os(4),phi1Os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NOsB_eu_a,NOsB_eu_b] ...
    = solvemodel(bhat1Os(1),[bhat1Os(2);0],bhat1Os(3),bhat1Os(4),phi1Os,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NOsC_eu_a,NOsC_eu_b] ...
    = solvemodel(bhat1Os(1),[bhat1Os(2);0],bhat1Os(3),bhat1Os(4),phi1Os,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), Without Omega(o), With Sigma(S)
    [~,~,~,~,~,NoSA_eu_a,NoSA_eu_b] ...
    = solvemodel(bhat1oS(1),[bhat1oS(2);0],bhat1oS(3),bhat1oS(4),phi1oS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NoSB_eu_a,NoSB_eu_b] ...
    = solvemodel(bhat1oS(1),[bhat1oS(2);0],bhat1oS(3),bhat1oS(4),phi1oS,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NoSC_eu_a,NoSC_eu_b] ...
    = solvemodel(bhat1oS(1),[bhat1oS(2);0],bhat1oS(3),bhat1oS(4),phi1oS,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), Without Omega(o), Without Sigma(s)
    [~,~,~,~,~,NosA_eu_a,NosA_eu_b] ...
    = solvemodel(bhat1os(1),[bhat1os(2);0],bhat1os(3),bhat1os(4),phi1os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NosB_eu_a,NosB_eu_b] ...
    = solvemodel(bhat1os(1),[bhat1os(2);0],bhat1os(3),bhat1os(4),phi1os,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    [~,~,~,~,~,NosC_eu_a,NosC_eu_b] ...
    = solvemodel(bhat1os(1),[bhat1os(2);0],bhat1os(3),bhat1os(4),phi1os,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

%% Sophisticated
sophisticated=1;
 
    % Sophisticated(S), With Omega(O), With Sigma(S)
    [~,~,~,~,~,SOSA_eu_a,SOSA_eu_b] ...
    = solvemodel(bhat2OS(1),[bhat2OS(2);0],bhat2OS(3),bhat2OS(4),phi2OS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SOSB_eu_a,SOSB_eu_b] ...
    = solvemodel(bhat2OS(1),[bhat2OS(2);0],bhat2OS(3),bhat2OS(4),phi2OS,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SOSC_eu_a,SOSC_eu_b] ...
    = solvemodel(bhat2OS(1),[bhat2OS(2);0],bhat2OS(3),bhat2OS(4),phi2OS,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    % Sophisticated(S), With Omega(O), Without Sigma(s)
    [~,~,~,~,~,SOsA_eu_a,SOsA_eu_b] ...
    = solvemodel(bhat2Os(1),[bhat2Os(2);0],bhat2Os(3),bhat2Os(4),phi2Os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);
 
    [~,~,~,~,~,SOsB_eu_a,SOsB_eu_b] ...
    = solvemodel(bhat2Os(1),[bhat2Os(2);0],bhat2Os(3),bhat2Os(4),phi2Os,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SOsC_eu_a,SOsC_eu_b] ...
    = solvemodel(bhat2Os(1),[bhat2Os(2);0],bhat2Os(3),bhat2Os(4),phi2Os,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    % Sophisticated(S), Without Omega(o), With Sigma(S)
    [~,~,~,~,~,SoSA_eu_a,SoSA_eu_b] ...
    = solvemodel(bhat2oS(1),[bhat2oS(2);0],bhat2oS(3),bhat2oS(4),phi2oS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SoSB_eu_a,SoSB_eu_b] ...
    = solvemodel(bhat2oS(1),[bhat2oS(2);0],bhat2oS(3),bhat2oS(4),phi2oS,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SoSC_eu_a,SoSC_eu_b] ...
    = solvemodel(bhat2oS(1),[bhat2oS(2);0],bhat2oS(3),bhat2oS(4),phi2oS,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    % Sophisticated(S), Without Omega(o), Without Sigma(s)
    [~,~,~,~,~,SosA_eu_a,SosA_eu_b] ...
    = solvemodel(bhat2os(1),[bhat2os(2);0],bhat2os(3),bhat2os(4),phi2os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SosB_eu_a,SosB_eu_b] ...
    = solvemodel(bhat2os(1),[bhat2os(2);0],bhat2os(3),bhat2os(4),phi2os,a1_,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

    [~,~,~,~,~,SosC_eu_a,SosC_eu_b] ...
    = solvemodel(bhat2os(1),[bhat2os(2);0],bhat2os(3),bhat2os(4),phi2os,a1,a2_,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

%% find the averages for each question

%typeq_extended= 8127*100x1
typeq_extended = reshape( repmat(typeq,1,numberdraws),[],1);

%%
%% NOS_eu_a_extended= 8127*100x1
NOSA_eu_a_extended= reshape(NOSA_eu_a,[],1);
NOSA_eu_b_extended= reshape(NOSA_eu_b,[],1);
NOsA_eu_a_extended= reshape(NOsA_eu_a,[],1);
NOsA_eu_b_extended= reshape(NOsA_eu_b,[],1);
NoSA_eu_a_extended= reshape(NoSA_eu_a,[],1);
NoSA_eu_b_extended= reshape(NoSA_eu_b,[],1);
NosA_eu_a_extended= reshape(NosA_eu_a,[],1);
NosA_eu_b_extended= reshape(NosA_eu_b,[],1);
SOSA_eu_a_extended= reshape(SOSA_eu_a,[],1);
SOSA_eu_b_extended= reshape(SOSA_eu_b,[],1);
SOsA_eu_a_extended= reshape(SOsA_eu_a,[],1);
SOsA_eu_b_extended= reshape(SOsA_eu_b,[],1);
SoSA_eu_a_extended= reshape(SoSA_eu_a,[],1);
SoSA_eu_b_extended= reshape(SoSA_eu_b,[],1);
SosA_eu_a_extended= reshape(SosA_eu_a,[],1);
SosA_eu_b_extended= reshape(SosA_eu_b,[],1);

%temp_NOS_eu_a_extended= 2709X100
temp_NOSA_eu_a_extended=reshape(NOSA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOSA_eu_b_extended=reshape(NOSA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsA_eu_a_extended=reshape(NOsA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsA_eu_b_extended=reshape(NOsA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSA_eu_a_extended=reshape(NoSA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSA_eu_b_extended=reshape(NoSA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosA_eu_a_extended=reshape(NosA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosA_eu_b_extended=reshape(NosA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSA_eu_a_extended=reshape(SOSA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSA_eu_b_extended=reshape(SOSA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsA_eu_a_extended=reshape(SOsA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsA_eu_b_extended=reshape(SOsA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSA_eu_a_extended=reshape(SoSA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSA_eu_b_extended=reshape(SoSA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosA_eu_a_extended=reshape(SosA_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosA_eu_b_extended=reshape(SosA_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);

%plot_NOS_eu_a = 21 (means taken over draws and then subjects)
plot_NOSA_eu_a = nanmean(nanmean(reshape(temp_NOSA_eu_a_extended,21,[],numberdraws),3),2);
plot_NOSA_eu_b = nanmean(nanmean(reshape(temp_NOSA_eu_b_extended,21,[],numberdraws),3),2);
plot_NOsA_eu_a = nanmean(nanmean(reshape(temp_NOsA_eu_a_extended,21,[],numberdraws),3),2);
plot_NOsA_eu_b = nanmean(nanmean(reshape(temp_NOsA_eu_b_extended,21,[],numberdraws),3),2);
plot_NoSA_eu_a = nanmean(nanmean(reshape(temp_NoSA_eu_a_extended,21,[],numberdraws),3),2);
plot_NoSA_eu_b = nanmean(nanmean(reshape(temp_NoSA_eu_b_extended,21,[],numberdraws),3),2);
plot_NosA_eu_a = nanmean(nanmean(reshape(temp_NosA_eu_a_extended,21,[],numberdraws),3),2);
plot_NosA_eu_b = nanmean(nanmean(reshape(temp_NosA_eu_b_extended,21,[],numberdraws),3),2);

plot_SOSA_eu_a = nanmean(nanmean(reshape(temp_SOSA_eu_a_extended,21,[],numberdraws),3),2);
plot_SOSA_eu_b = nanmean(nanmean(reshape(temp_SOSA_eu_b_extended,21,[],numberdraws),3),2);
plot_SOsA_eu_a = nanmean(nanmean(reshape(temp_SOsA_eu_a_extended,21,[],numberdraws),3),2);
plot_SOsA_eu_b = nanmean(nanmean(reshape(temp_SOsA_eu_b_extended,21,[],numberdraws),3),2);
plot_SoSA_eu_a = nanmean(nanmean(reshape(temp_SoSA_eu_a_extended,21,[],numberdraws),3),2);
plot_SoSA_eu_b = nanmean(nanmean(reshape(temp_SoSA_eu_b_extended,21,[],numberdraws),3),2);
plot_SosA_eu_a = nanmean(nanmean(reshape(temp_SosA_eu_a_extended,21,[],numberdraws),3),2);
plot_SosA_eu_b = nanmean(nanmean(reshape(temp_SosA_eu_b_extended,21,[],numberdraws),3),2);

% discount by beta delta
delta1=mean(max(0,1-par_draws1(:,2).^2),1);
beta1=mean(max(0,1-par_draws1(:,1).^2),1);
delta2=mean(max(0,1-par_draws2(:,2).^2),1);
beta2=mean(max(0,1-par_draws2(:,1).^2),1);
plot_NOSA_eu_a = plot_NOSA_eu_a * delta1 * beta1; 
plot_NOSA_eu_b = plot_NOSA_eu_b * delta1 * beta1; 
plot_NOsA_eu_a = plot_NOsA_eu_a * delta1 * beta1; 
plot_NOsA_eu_b = plot_NOsA_eu_b * delta1 * beta1; 
plot_NoSA_eu_a = plot_NoSA_eu_a * delta1 * beta1; 
plot_NoSA_eu_b = plot_NoSA_eu_b * delta1 * beta1; 
plot_NosA_eu_a = plot_NosA_eu_a * delta1 * beta1; 
plot_NosA_eu_b = plot_NosA_eu_b * delta1 * beta1; 

plot_SOSA_eu_a = plot_SOSA_eu_a * delta2 * beta2; 
plot_SOSA_eu_b = plot_SOSA_eu_b * delta2 * beta2; 
plot_SOsA_eu_a = plot_SOsA_eu_a * delta2 * beta2; 
plot_SOsA_eu_b = plot_SOsA_eu_b * delta2 * beta2; 
plot_SoSA_eu_a = plot_SoSA_eu_a * delta2 * beta2; 
plot_SoSA_eu_b = plot_SoSA_eu_b * delta2 * beta2;
plot_SosA_eu_a = plot_SosA_eu_a * delta2 * beta2;
plot_SosA_eu_b = plot_SosA_eu_b * delta2 * beta2;

%%
%% NOS_eu_a_extended= 8127*100x1
NOSB_eu_a_extended= reshape(NOSB_eu_a,[],1);
NOSB_eu_b_extended= reshape(NOSB_eu_b,[],1);
NOsB_eu_a_extended= reshape(NOsB_eu_a,[],1);
NOsB_eu_b_extended= reshape(NOsB_eu_b,[],1);
NoSB_eu_a_extended= reshape(NoSB_eu_a,[],1);
NoSB_eu_b_extended= reshape(NoSB_eu_b,[],1);
NosB_eu_a_extended= reshape(NosB_eu_a,[],1);
NosB_eu_b_extended= reshape(NosB_eu_b,[],1);
SOSB_eu_a_extended= reshape(SOSB_eu_a,[],1);
SOSB_eu_b_extended= reshape(SOSB_eu_b,[],1);
SOsB_eu_a_extended= reshape(SOsB_eu_a,[],1);
SOsB_eu_b_extended= reshape(SOsB_eu_b,[],1);
SoSB_eu_a_extended= reshape(SoSB_eu_a,[],1);
SoSB_eu_b_extended= reshape(SoSB_eu_b,[],1);
SosB_eu_a_extended= reshape(SosB_eu_a,[],1);
SosB_eu_b_extended= reshape(SosB_eu_b,[],1);

%temp_NOS_eu_a_extended= 2709X100
temp_NOSB_eu_a_extended=reshape(NOSB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOSB_eu_b_extended=reshape(NOSB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsB_eu_a_extended=reshape(NOsB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsB_eu_b_extended=reshape(NOsB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSB_eu_a_extended=reshape(NoSB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSB_eu_b_extended=reshape(NoSB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosB_eu_a_extended=reshape(NosB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosB_eu_b_extended=reshape(NosB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSB_eu_a_extended=reshape(SOSB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSB_eu_b_extended=reshape(SOSB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsB_eu_a_extended=reshape(SOsB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsB_eu_b_extended=reshape(SOsB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSB_eu_a_extended=reshape(SoSB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSB_eu_b_extended=reshape(SoSB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosB_eu_a_extended=reshape(SosB_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosB_eu_b_extended=reshape(SosB_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);

%plot_NOS_eu_a = 21 (means taken over draws and then subjects)
plot_NOSB_eu_a = nanmean(nanmean(reshape(temp_NOSB_eu_a_extended,21,[],numberdraws),3),2);
plot_NOSB_eu_b = nanmean(nanmean(reshape(temp_NOSB_eu_b_extended,21,[],numberdraws),3),2);
plot_NOsB_eu_a = nanmean(nanmean(reshape(temp_NOsB_eu_a_extended,21,[],numberdraws),3),2);
plot_NOsB_eu_b = nanmean(nanmean(reshape(temp_NOsB_eu_b_extended,21,[],numberdraws),3),2);
plot_NoSB_eu_a = nanmean(nanmean(reshape(temp_NoSB_eu_a_extended,21,[],numberdraws),3),2);
plot_NoSB_eu_b = nanmean(nanmean(reshape(temp_NoSB_eu_b_extended,21,[],numberdraws),3),2);
plot_NosB_eu_a = nanmean(nanmean(reshape(temp_NosB_eu_a_extended,21,[],numberdraws),3),2);
plot_NosB_eu_b = nanmean(nanmean(reshape(temp_NosB_eu_b_extended,21,[],numberdraws),3),2);

plot_SOSB_eu_a = nanmean(nanmean(reshape(temp_SOSB_eu_a_extended,21,[],numberdraws),3),2);
plot_SOSB_eu_b = nanmean(nanmean(reshape(temp_SOSB_eu_b_extended,21,[],numberdraws),3),2);
plot_SOsB_eu_a = nanmean(nanmean(reshape(temp_SOsB_eu_a_extended,21,[],numberdraws),3),2);
plot_SOsB_eu_b = nanmean(nanmean(reshape(temp_SOsB_eu_b_extended,21,[],numberdraws),3),2);
plot_SoSB_eu_a = nanmean(nanmean(reshape(temp_SoSB_eu_a_extended,21,[],numberdraws),3),2);
plot_SoSB_eu_b = nanmean(nanmean(reshape(temp_SoSB_eu_b_extended,21,[],numberdraws),3),2);
plot_SosB_eu_a = nanmean(nanmean(reshape(temp_SosB_eu_a_extended,21,[],numberdraws),3),2);
plot_SosB_eu_b = nanmean(nanmean(reshape(temp_SosB_eu_b_extended,21,[],numberdraws),3),2);

% discount by beta delta
delta1=mean(max(0,1-par_draws1(:,2).^2),1);
beta1=mean(max(0,1-par_draws1(:,1).^2),1);
delta2=mean(max(0,1-par_draws2(:,2).^2),1);
beta2=mean(max(0,1-par_draws2(:,1).^2),1);
plot_NOSB_eu_a = plot_NOSB_eu_a * delta1 * beta1; 
plot_NOSB_eu_b = plot_NOSB_eu_b * delta1 * beta1; 
plot_NOsB_eu_a = plot_NOsB_eu_a * delta1 * beta1; 
plot_NOsB_eu_b = plot_NOsB_eu_b * delta1 * beta1; 
plot_NoSB_eu_a = plot_NoSB_eu_a * delta1 * beta1; 
plot_NoSB_eu_b = plot_NoSB_eu_b * delta1 * beta1; 
plot_NosB_eu_a = plot_NosB_eu_a * delta1 * beta1; 
plot_NosB_eu_b = plot_NosB_eu_b * delta1 * beta1; 

plot_SOSB_eu_a = plot_SOSB_eu_a * delta2 * beta2; 
plot_SOSB_eu_b = plot_SOSB_eu_b * delta2 * beta2; 
plot_SOsB_eu_a = plot_SOsB_eu_a * delta2 * beta2; 
plot_SOsB_eu_b = plot_SOsB_eu_b * delta2 * beta2; 
plot_SoSB_eu_a = plot_SoSB_eu_a * delta2 * beta2; 
plot_SoSB_eu_b = plot_SoSB_eu_b * delta2 * beta2;
plot_SosB_eu_a = plot_SosB_eu_a * delta2 * beta2;
plot_SosB_eu_b = plot_SosB_eu_b * delta2 * beta2;

%%
%% NOS_eu_a_extended= 8127*100x1
NOSC_eu_a_extended= reshape(NOSC_eu_a,[],1);
NOSC_eu_b_extended= reshape(NOSC_eu_b,[],1);
NOsC_eu_a_extended= reshape(NOsC_eu_a,[],1);
NOsC_eu_b_extended= reshape(NOsC_eu_b,[],1);
NoSC_eu_a_extended= reshape(NoSC_eu_a,[],1);
NoSC_eu_b_extended= reshape(NoSC_eu_b,[],1);
NosC_eu_a_extended= reshape(NosC_eu_a,[],1);
NosC_eu_b_extended= reshape(NosC_eu_b,[],1);
SOSC_eu_a_extended= reshape(SOSC_eu_a,[],1);
SOSC_eu_b_extended= reshape(SOSC_eu_b,[],1);
SOsC_eu_a_extended= reshape(SOsC_eu_a,[],1);
SOsC_eu_b_extended= reshape(SOsC_eu_b,[],1);
SoSC_eu_a_extended= reshape(SoSC_eu_a,[],1);
SoSC_eu_b_extended= reshape(SoSC_eu_b,[],1);
SosC_eu_a_extended= reshape(SosC_eu_a,[],1);
SosC_eu_b_extended= reshape(SosC_eu_b,[],1);

%temp_NOS_eu_a_extended= 2709X100
temp_NOSC_eu_a_extended=reshape(NOSC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOSC_eu_b_extended=reshape(NOSC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsC_eu_a_extended=reshape(NOsC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOsC_eu_b_extended=reshape(NOsC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSC_eu_a_extended=reshape(NoSC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoSC_eu_b_extended=reshape(NoSC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosC_eu_a_extended=reshape(NosC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NosC_eu_b_extended=reshape(NosC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSC_eu_a_extended=reshape(SOSC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOSC_eu_b_extended=reshape(SOSC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsC_eu_a_extended=reshape(SOsC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOsC_eu_b_extended=reshape(SOsC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSC_eu_a_extended=reshape(SoSC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoSC_eu_b_extended=reshape(SoSC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosC_eu_a_extended=reshape(SosC_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SosC_eu_b_extended=reshape(SosC_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);

%plot_NOS_eu_a = 21 (means taken over draws and then subjects)
plot_NOSC_eu_a = nanmean(nanmean(reshape(temp_NOSC_eu_a_extended,21,[],numberdraws),3),2);
plot_NOSC_eu_b = nanmean(nanmean(reshape(temp_NOSC_eu_b_extended,21,[],numberdraws),3),2);
plot_NOsC_eu_a = nanmean(nanmean(reshape(temp_NOsC_eu_a_extended,21,[],numberdraws),3),2);
plot_NOsC_eu_b = nanmean(nanmean(reshape(temp_NOsC_eu_b_extended,21,[],numberdraws),3),2);
plot_NoSC_eu_a = nanmean(nanmean(reshape(temp_NoSC_eu_a_extended,21,[],numberdraws),3),2);
plot_NoSC_eu_b = nanmean(nanmean(reshape(temp_NoSC_eu_b_extended,21,[],numberdraws),3),2);
plot_NosC_eu_a = nanmean(nanmean(reshape(temp_NosC_eu_a_extended,21,[],numberdraws),3),2);
plot_NosC_eu_b = nanmean(nanmean(reshape(temp_NosC_eu_b_extended,21,[],numberdraws),3),2);

plot_SOSC_eu_a = nanmean(nanmean(reshape(temp_SOSC_eu_a_extended,21,[],numberdraws),3),2);
plot_SOSC_eu_b = nanmean(nanmean(reshape(temp_SOSC_eu_b_extended,21,[],numberdraws),3),2);
plot_SOsC_eu_a = nanmean(nanmean(reshape(temp_SOsC_eu_a_extended,21,[],numberdraws),3),2);
plot_SOsC_eu_b = nanmean(nanmean(reshape(temp_SOsC_eu_b_extended,21,[],numberdraws),3),2);
plot_SoSC_eu_a = nanmean(nanmean(reshape(temp_SoSC_eu_a_extended,21,[],numberdraws),3),2);
plot_SoSC_eu_b = nanmean(nanmean(reshape(temp_SoSC_eu_b_extended,21,[],numberdraws),3),2);
plot_SosC_eu_a = nanmean(nanmean(reshape(temp_SosC_eu_a_extended,21,[],numberdraws),3),2);
plot_SosC_eu_b = nanmean(nanmean(reshape(temp_SosC_eu_b_extended,21,[],numberdraws),3),2);

% discount by beta delta
delta1=mean(max(0,1-par_draws1(:,2).^2),1);
beta1=mean(max(0,1-par_draws1(:,1).^2),1);
delta2=mean(max(0,1-par_draws2(:,2).^2),1);
beta2=mean(max(0,1-par_draws2(:,1).^2),1);
plot_NOSC_eu_a = plot_NOSC_eu_a * delta1 * beta1; 
plot_NOSC_eu_b = plot_NOSC_eu_b * delta1 * beta1; 
plot_NOsC_eu_a = plot_NOsC_eu_a * delta1 * beta1; 
plot_NOsC_eu_b = plot_NOsC_eu_b * delta1 * beta1; 
plot_NoSC_eu_a = plot_NoSC_eu_a * delta1 * beta1; 
plot_NoSC_eu_b = plot_NoSC_eu_b * delta1 * beta1; 
plot_NosC_eu_a = plot_NosC_eu_a * delta1 * beta1; 
plot_NosC_eu_b = plot_NosC_eu_b * delta1 * beta1; 

plot_SOSC_eu_a = plot_SOSC_eu_a * delta2 * beta2; 
plot_SOSC_eu_b = plot_SOSC_eu_b * delta2 * beta2; 
plot_SOsC_eu_a = plot_SOsC_eu_a * delta2 * beta2; 
plot_SOsC_eu_b = plot_SOsC_eu_b * delta2 * beta2; 
plot_SoSC_eu_a = plot_SoSC_eu_a * delta2 * beta2; 
plot_SoSC_eu_b = plot_SoSC_eu_b * delta2 * beta2;
plot_SosC_eu_a = plot_SosC_eu_a * delta2 * beta2;
plot_SosC_eu_b = plot_SosC_eu_b * delta2 * beta2;

%% payment amounts
payments = [ 4  1  6      8.24    NaN      7.37 ;
    4  2  6      8.88    NaN      7.94 ;
    4  3  6      10.14   NaN      9.07 ;
    4  4  6      11.41   NaN      10.21 ;
    4  5  6      12.68   NaN      11.34 ;
    4  6  6      13.95   NaN      12.47 ;
    4  7  6      15.21   NaN      13.61 ;
    5  8  6      8.24   6.5  NaN   ;
    5  9  6      8.88   7    NaN   ;
    5  10  6      10.14  8   NaN    ;
    5  11  6      11.41  9   NaN    ;
    5  12  6      12.68  10  NaN    ;
    5  13  6      13.95  11  NaN    ;
    5  14  6      15.21  12  NaN    ;
    6  15  6      8.24   6.5    7.37 ;
    6  16  6      8.88   7      7.94 ;
    6  17  6      10.14  8      9.07 ;
    6  18  6      11.41  9      10.21 ;
    6  19  6      12.68  10     11.34 ;
    6  20  6      13.95  11     12.47 ;
    6  21  6      15.21  12     13.61 ];

%% tables for 4 categories NO , No, SO, So

fprintf('%% Table 2, Naive with omega \n')

input.data = [payments plot_NOSA_eu_a plot_NOsA_eu_a plot_NOSB_eu_a plot_NOsB_eu_a plot_NOSC_eu_a plot_NOsC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A' '$\sigma = 0$, A' '$\sigma \neq 0$, B' '$\sigma = 0$, B' '$\sigma \neq 0$, C' '$\sigma = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\omega \neq 0$';
% LaTex table label:
input.tableLabel = 'naivewithomega2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%No
fprintf('%% Table 2, Naive without omega \n')
input.data = [payments plot_NoSA_eu_a plot_NosA_eu_a plot_NoSB_eu_a plot_NosB_eu_a plot_NoSC_eu_a plot_NosC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A' '$\sigma = 0$, A' '$\sigma \neq 0$, B' '$\sigma = 0$, B' '$\sigma \neq 0$, C' '$\sigma = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\omega = 0$';
% LaTex table label:
input.tableLabel = 'naivewithoutomega2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%SO
fprintf('%% Table 2, Sophisticated with omega \n')
input.data = [payments plot_SOSA_eu_a plot_SOsA_eu_a plot_SOSB_eu_a plot_SOsB_eu_a plot_SOSC_eu_a plot_SOsC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A' '$\sigma = 0$, A' '$\sigma \neq 0$, B' '$\sigma = 0$, B' '$\sigma \neq 0$, C' '$\sigma = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\omega \neq 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithomega2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%So
fprintf('%% Table 2, Sophisticated without omega \n')
input.data = [payments plot_SoSA_eu_a plot_SosA_eu_a plot_SoSB_eu_a plot_SosB_eu_a plot_SoSC_eu_a plot_SosC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A' '$\sigma = 0$, A' '$\sigma \neq 0$, B' '$\sigma = 0$, B' '$\sigma \neq 0$, C' '$\sigma = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\omega = 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithoutomega2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%% tables for 4 categories NS , Ns, SS, Ss

% NS

fprintf('%% Table 2, Naive with sigma \n')

input.data = [payments plot_NOSA_eu_a plot_NoSA_eu_a plot_NOSB_eu_a plot_NoSB_eu_a plot_NOSC_eu_a plot_NoSC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A' '$\omega = 0$, A' '$\omega \neq 0$, B' '$\omega = 0$, B' '$\omega \neq 0$, C' '$\omega = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\sigma \neq 0$';
% LaTex table label:
input.tableLabel = 'naivewithsigma2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% Ns

fprintf('%% Table 2, Naive without sigma \n')

input.data = [payments plot_NOsA_eu_a plot_NosA_eu_a plot_NOsB_eu_a plot_NosB_eu_a plot_NOsC_eu_a plot_NosC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A' '$\omega = 0$, A' '$\omega \neq 0$, B' '$\omega = 0$, B' '$\omega \neq 0$, C' '$\omega = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\sigma = 0$';
% LaTex table label:
input.tableLabel = 'naivewithoutsigma2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% SS

fprintf('%% Table 2, Sophisticated with sigma \n')

input.data = [payments plot_SOSA_eu_a plot_SoSA_eu_a plot_SOSB_eu_a plot_SoSB_eu_a plot_SOSC_eu_a plot_SoSC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A' '$\omega = 0$, A' '$\omega \neq 0$, B' '$\omega = 0$, B' '$\omega \neq 0$, C' '$\omega = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\sigma \neq 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithsigma2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% Ss

fprintf('%% Table 2, Sophisticated without sigma \n')

input.data = [payments plot_SOsA_eu_a plot_SosA_eu_a plot_SOsB_eu_a plot_SosB_eu_a plot_SOsC_eu_a plot_SosC_eu_a];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A' '$\omega = 0$, A' '$\omega \neq 0$, B' '$\omega = 0$, B' '$\omega \neq 0$, C' '$\omega = 0$, C'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',10}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\sigma = 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithoutsigma2';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

fprintf('%% END OF COPY AND PASTE \n')

%% figures for 4 categories NO , No, SO, So

%NO
close all
subplot(1,2,1) 
bar(1:21,[plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]))]);
hleg1 = legend('A \sigma \neq 0','B \sigma \neq 0','C \sigma \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]))]);
hleg1 = legend('A \sigma = 0','B \sigma = 0','C \sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega \neq 0')  

print('-dpdf','../output/matlab_table2_NbigO2.pdf')

%No
close all
subplot(1,2,1) 
bar(1:21,[plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
hleg1 = legend('A \sigma \neq 0','B \sigma \neq 0','C \sigma \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega = 0')  

subplot(1,2,2) 
bar(1:21,[plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
hleg1 = legend('A \sigma = 0','B \sigma = 0','C \sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega = 0')

print('-dpdf','../output/matlab_table2_Nsmallo2.pdf')

%SO
close all
subplot(1,2,1) 
bar(1:21,[plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]))]);
hleg1 = legend('A \sigma \neq 0','B \sigma \neq 0','C \sigma \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]))]);
hleg1 = legend('A \sigma = 0','B \sigma = 0','C \sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega \neq 0')  

print('-dpdf','../output/matlab_table2_SbigO2.pdf')

%So
close all
subplot(1,2,1) 
bar(1:21,[plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
hleg1 = legend('A \sigma \neq 0','B \sigma \neq 0','C \sigma \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega = 0')  

subplot(1,2,2) 
bar(1:21,[plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
hleg1 = legend('A \sigma = 0','B \sigma = 0','C \sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega = 0')  

print('-dpdf','../output/matlab_table2_Ssmallo.pdf')

%% figures for 4 categories NS , Ns, SS, Ss

% NS
close all
subplot(1,2,1) 
bar(1:21,[plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]))]);
hleg1 = legend('A \omega \neq 0','B \omega \neq 0','C \omega \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOSA_eu_a plot_NOSB_eu_a plot_NOSC_eu_a plot_NoSA_eu_a plot_NoSB_eu_a plot_NoSC_eu_a]))]);
hleg1 = legend('A \omega = 0','B \omega = 0','C \omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma \neq 0')  

print('-dpdf','../output/matlab_table2_NbigS2.pdf')

% Ns
close all
subplot(1,2,1) 
bar(1:21,[plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
hleg1 = legend('A \omega \neq 0','B \omega \neq 0','C \omega \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma = 0')  

subplot(1,2,2) 
bar(1:21,[plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOsA_eu_a plot_NOsB_eu_a plot_NOsC_eu_a plot_NosA_eu_a plot_NosB_eu_a plot_NosC_eu_a]))]);
hleg1 = legend('A \omega = 0','B \omega = 0','C \omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma = 0')  

print('-dpdf','../output/matlab_table2_NsmallS.pdf')

% SS
close all
subplot(1,2,1) 
bar(1:21,[plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]))]);
hleg1 = legend('A \omega \neq 0','B \omega \neq 0','C \omega \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOSA_eu_a plot_SOSB_eu_a plot_SOSC_eu_a plot_SoSA_eu_a plot_SoSB_eu_a plot_SoSC_eu_a]))]);
hleg1 = legend('A \omega = 0','B \omega = 0','C \omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma \neq 0')  

print('-dpdf','../output/matlab_table2_SbigS2.pdf')

% Ss
close all
subplot(1,2,1) 
bar(1:21,[plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
hleg1 = legend('A \omega \neq 0','B \omega \neq 0','C \omega \neq 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma = 0')  

subplot(1,2,2) 
bar(1:21,[plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]);
barmap=[255/256 255/256 0/256; 0/256 0/256 255/256; 255/256 0/256 0/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOsA_eu_a plot_SOsB_eu_a plot_SOsC_eu_a plot_SosA_eu_a plot_SosB_eu_a plot_SosC_eu_a]))]);
hleg1 = legend('A \omega = 0','B \omega = 0','C \omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma = 0')  

print('-dpdf','../output/matlab_table2_Ssmalls2.pdf')

%% example 05052015 can be commented out

results=[ 
mean(plot_NOSC_eu_a,1) mean(plot_NOSB_eu_a,1) mean(plot_NOSA_eu_a,1); ...
mean(plot_NOsC_eu_a,1) mean(plot_NOsB_eu_a,1) mean(plot_NOsA_eu_a,1); ...
mean(plot_NoSC_eu_a,1) mean(plot_NoSB_eu_a,1) mean(plot_NoSA_eu_a,1); ...
mean(plot_NosC_eu_a,1) mean(plot_NoSB_eu_a,1) mean(plot_NoSA_eu_a,1); ...
mean(plot_SOSC_eu_a,1) mean(plot_SOSB_eu_a,1) mean(plot_SOSA_eu_a,1); ...
mean(plot_SOsC_eu_a,1) mean(plot_SOsB_eu_a,1) mean(plot_SOsA_eu_a,1); ...
mean(plot_SoSC_eu_a,1) mean(plot_SoSB_eu_a,1) mean(plot_SoSA_eu_a,1); ...
mean(plot_SosC_eu_a,1) mean(plot_SosB_eu_a,1) mean(plot_SosA_eu_a,1); ...
]

%results=[NaN(size(results,1),size(results,2)) , results];

input.data = results;
input.tableColLabels = {  '$m_{ije}=10$' '$m_{ij\ell}=10/\bar{\delta}$' '$m_{ije}=10$, $m_{ij\ell}=10/\bar{\delta}$'};
% input.tableRowLabels = {'label1' 'label2' };
input.tableRowLabels = {'Naive, $\omega \neq 0$, $\sigma \neq 0$' ...
'Naive, $\omega \neq 0$, $\sigma = 0$' ...
'Naive, $\omega = 0$, $\sigma \neq 0$' ...
'Naive, $\omega = 0$, $\sigma = 0$' ...
'Sophisticated, $\omega \neq 0$, $\sigma \neq 0$' ...
'Sophisticated, $\omega \neq 0$, $\sigma = 0$' ...
'Sophisticated, $\omega = 0$, $\sigma \neq 0$' ...
'Sophisticated, $\omega = 0$, $\sigma = 0$' ...
    };
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.2f',3}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'The average value of receiving 10 Euros two months later from the perspective of self-0';
% LaTex table label:
input.tableLabel = 'exercise';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

