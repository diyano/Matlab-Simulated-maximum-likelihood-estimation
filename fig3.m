%% figure(3)
close all

%21 question fit
[~,~,~,~,N_stage1_prob_a,S_stage1_prob_a,~,~] = objfun(bhat_fmincon);
num_dim=size(startvalues,1);

typeq_extended = reshape( repmat(typeq,1,numberdraws),[],1);
temp_N_stage1_prob_a = reshape(N_stage1_prob_a,[],1);
temp_S_stage1_prob_a = reshape(S_stage1_prob_a,[],1);
temp_N_stage1_prob_a = reshape(temp_N_stage1_prob_a(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_S_stage1_prob_a = reshape(temp_S_stage1_prob_a(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);

fit_model_raw=reshape(temp_N_stage1_prob_a,21,[],numberdraws) .* (1-bhat_fmincon(num_dim)) + reshape(temp_S_stage1_prob_a,21,[],numberdraws) .* bhat_fmincon(num_dim);
fit_model = mean(mean(fit_model_raw,3),2);

temp_choose_a= choose_a(typeq==4|typeq==5|typeq==6);
fit_data_raw = reshape(temp_choose_a,21,[]);
fit_data = mean(fit_data_raw,2);

%% figure(3) naive

fit_model_raw=reshape(temp_N_stage1_prob_a,21,[],numberdraws); 
fit_model_naive = mean(mean(fit_model_raw,3),2);


%% figure(3) sophisticated

fit_model_raw=reshape(temp_S_stage1_prob_a,21,[],numberdraws);
fit_model_sophisticated = mean(mean(fit_model_raw,3),2);

% figure exporting options
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

subplot(1,1,1) 
bar(1:21,[1-fit_data 1-fit_model])
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([1-fit_data 1-fit_model_naive 1-fit_model_sophisticated]))])
set(gca,'XTick',1:21)
hleg1 = legend('Data','Model');
xlabel('Question','FontSize',fontsize)
ylabel('Probability','FontSize',fontsize)
set(hleg1,'Location','NorthWest')

print('-dpdf','../output/matlab_figure3.pdf')

subplot(1,1,1)
bar(1:21,[1-fit_data 1-fit_model_naive 1-fit_model_sophisticated])
barmap=[207/256 181/256 59/256; 173/256 216/256 230/256; 0 0 139/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([1-fit_data 1-fit_model_naive 1-fit_model_sophisticated]))])
set(gca,'XTick',1:21)
hleg1 = legend('Data','Model, Naive', 'Model, Sophisticated');
xlabel('Question','FontSize',fontsize)
ylabel('Probability','FontSize',fontsize)
set(hleg1,'Location','NorthWest')

print('-dpdf','../output/matlab_figure3ns.pdf')