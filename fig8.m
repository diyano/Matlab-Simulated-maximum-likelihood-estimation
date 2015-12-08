% figure(8)
close all

delta1=max(0,1-par_draws1(:,2).^2);
beta1=min(bhat1(1) + 4.*sqrt(rc_vcov1(1)),max(bhat1(1) - 4.*sqrt(rc_vcov1(1)),par_draws1(:,1)));
delta2=max(0,1-par_draws2(:,2).^2);
beta2=min(bhat2(1) + 4.*sqrt(rc_vcov2(1)),max(bhat2(1) - 4.*sqrt(rc_vcov2(1)),par_draws2(:,1)));

num_intervals=20;
step_delta1=(max(delta1) - min(delta1))/num_intervals;
step_beta1=(max(beta1) - min(beta1))/num_intervals;
step_delta2=(max(delta2) - min(delta2))/num_intervals;
step_beta2=(max(beta2) - min(beta2))/num_intervals;

delta1_intervals=(min(delta1):step_delta1:max(delta1));
beta1_intervals=(min(beta1):step_beta1:max(beta1));
delta2_intervals=(min(delta2):step_delta2:max(delta2));
beta2_intervals=(min(beta2):step_beta2:max(beta2));

delta1_=(min(delta1)+step_delta1/2:step_delta1:max(delta1)-step_delta1/2);
beta1_=(min(beta1)+step_beta1/2:step_beta1:max(beta1)-step_beta1/2);
delta2_=(min(delta2)+step_delta2/2:step_delta2:max(delta2)-step_delta2/2);
beta2_=(min(beta2)+step_beta2/2:step_beta2:max(beta2)-step_beta2/2);

%% naive
counts1=NaN(num_intervals,num_intervals);
area1=[];
delta1_plot=[];
beta1_plot=[];

for i=1:num_intervals    
    for m=1:num_intervals
    n=m+1;
    j=i+1;
    counts1(i,m)=size(delta1(delta1>=delta1_intervals(i) & delta1<delta1_intervals(j) & beta1 >= beta1_intervals(m) & beta1<beta1_intervals(n)),1); 
        if counts1(i,m)>0
        delta1_plot=[delta1_plot,delta1_(i)];
        beta1_plot=[beta1_plot,beta1_(m)];
        area1=[area1,counts1(i,m)];
        end    
    end
end

p = polyfit(delta1,beta1,1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
fitted1 = p(1) .* delta1_plot + p(2); % compute a new vector r that has matching datapoints in x

subplot(1,1,1)
scatter(delta1_plot,beta1_plot,area1)
hold on
plot(delta1_plot,fitted1);
hold off
xlabel('\delta_i','FontSize',fontsize)
ylabel('\beta_i','FontSize',fontsize)
%set(gca,'yticklabel',num2str(get(gca,'ytick')','%4.10f'))
%title('Naive')

print('-dpdf','../output/matlab_figure8n.pdf')

%% sophisticated

counts2=NaN(num_intervals,num_intervals);
area2=[];
delta2_plot=[];
beta2_plot=[];

for i=1:num_intervals    
    for m=1:num_intervals
    n=m+1;
    j=i+1;
    counts2(i,m)=size(delta2(delta2>=delta2_intervals(i) & delta2<delta2_intervals(j) & beta2 >= beta2_intervals(m) & beta2<beta2_intervals(n)),1); 
        if counts2(i,m)>0
        delta2_plot=[delta2_plot,delta2_(i)];
        beta2_plot=[beta2_plot,beta2_(m)];
        area2=[area2,counts2(i,m)];
        end    
    end
end

p = polyfit(delta2,beta2,1);   
fitted2 = p(1) .* delta2_plot + p(2); 

subplot(1,1,1)
scatter(delta2_plot,beta2_plot,area2)
hold on
plot(delta2_plot,fitted2);
hold off
xlabel('\delta_i','FontSize',fontsize)
ylabel('\beta_i','FontSize',fontsize)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%4.10f'))
%title('Naive')

print('-dpdf','../output/matlab_figure8s.pdf')