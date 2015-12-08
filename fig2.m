%% plot distribution random coefficients
close all
% figure(2) naive
% figure exporting options
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

subplot(2,2,1)
draws_beta = par_draws1(:,1);
%histogram(draws_beta,'Normalization','probability')
hist(draws_beta)
xlabel('\beta_i','FontSize',fontsize)
ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,2)
draws_delta = max(0,min(1-par_draws1(:,2).^2,0.9999));
% histogram(draws_delta,'Normalization','probability')
hist(draws_delta)
xlabel('\delta_i','FontSize',fontsize)
%ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,3)
draws_sigma = max(1e-1,min(par_draws1(:,3).^2,10));
%histogram(draws_sigma,'Normalization','probability')
hist(draws_sigma)
xlabel('\sigma_i','FontSize',fontsize)
ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,4)
draws_omega = par_draws1(:,4);
%histogram(draws_omega,'Normalization','probability')
hist(draws_omega)
xlabel('\omega_i','FontSize',fontsize)
%ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

print('-dpdf','../output/matlab_figure2n.pdf')

% figure(2) sophisticated
% figure exporting options
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

par_draws=mvnrnd( bhat2',rc_vcov2,1000);

subplot(2,2,1)
draws_beta = par_draws2(:,1);
%histogram(draws_beta,'Normalization','probability')
hist(draws_beta)
xlabel('\beta_i','FontSize',fontsize)
ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,2)
draws_delta = max(0,min(1-par_draws2(:,2).^2,0.9999));
%histogram(draws_delta,'Normalization','probability')
hist(draws_delta)
xlabel('\delta_i','FontSize',fontsize)
%ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,3)
draws_sigma = max(1e-1,min(par_draws2(:,3).^2,10));
%histogram(draws_sigma,'Normalization','probability')
hist(draws_sigma)
xlabel('\sigma_i','FontSize',fontsize)
ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

subplot(2,2,4)
draws_omega = par_draws2(:,4);
%histogram(draws_omega,'Normalization','probability')
hist(draws_omega)
xlabel('\omega_i','FontSize',fontsize)
%ylabel('Frequency','FontSize',fontsize)
set(gca, 'YTick', []);

print('-dpdf','../output/matlab_figure2s.pdf')
