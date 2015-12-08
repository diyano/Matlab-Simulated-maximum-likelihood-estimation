% figure(4)
close all

[~,~,~,~,N_stage1_prob_a,S_stage1_prob_a,stage2_prob_early_ifa,stage2_prob_early_ifb] = objfun(bhat_fmincon);
num_dim=size(startvalues,1);

% figure exporting options
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

% first mean is taken across draws and the second over questions
fit_model_a = mean(reshape(mean(stage2_prob_early_ifa(isnan(stage2_prob_early_ifa)==0),2),21,[]),2);
fit_model_b = mean(reshape(mean(stage2_prob_early_ifb(isnan(stage2_prob_early_ifb)==0),2),21,[]),2);

% choose_a : whether you choose a or not in the first stage (8127 x 1)  
% choose_a_stage2 : the answer given in second stage (1 late 0 early n/a no answer):

temp=reshape(choose_a.*choose_a_stage2,21,[]);
temp(logical(1-choose_a))=NaN;
temp2=temp;
temp3=1-isnan(temp);
temp2(isnan(temp2))=0;
fit_data_a = sum(temp2,2)./sum(temp3,2);

temp=reshape((1-choose_a).*choose_a_stage2,21,[]);
temp(logical(choose_a))=NaN;
temp2=temp;
temp3=1-isnan(temp);
temp2(isnan(temp2))=0;
fit_data_b = sum(temp2,2)./sum(temp3,2);
fit_data_b(1:7) = 0;
fit_data_b(8:14) = 1;

fit_model_a(isnan(fit_data_a))=NaN;
fit_model_b(isnan(fit_data_b))=NaN;

subplot(1,1,1)
bar(1:21,[1-fit_data_a 1-fit_model_a])
barmap=[207/256 181/256 59/256; 173/256 216/256 230/256; 0 0 139/256]; 
colormap(barmap);
axis([0 22 0 1])
set(gca,'XTick',1:21)
xlabel('Question','FontSize',fontsize)
ylabel('Probability','FontSize',fontsize)
print('-dpdf','../output/matlab_figure4a.pdf')

subplot(1,1,1)
bar(1:21,[1-fit_data_b 1-fit_model_b])
barmap=[207/256 181/256 59/256; 173/256 216/256 230/256; 0 0 139/256]; 
colormap(barmap);
axis([0 22 0 1])
set(gca,'XTick',1:21)
hleg1 = legend('data','model');
xlabel('Question','FontSize',fontsize)
ylabel('Probability','FontSize',fontsize)
print('-dpdf','../output/matlab_figure4b.pdf')