% random starting values or not. 1 for random starting values
random=1; 
numbermultistart = 12; % how many trials if random
% initialize matrices
bhat_multi = [];
nll_multi = [];
startvalues_multi= [];

if random==1 % if random starting values
for i=1:numbermultistart
    startvalues_random = mvnrnd(bhat_fmincon,diag(se_fmincon));
    startvalues_random(end)= min(0.9999,max(0.0001,startvalues_random(end)));
    [temp_bhat,temp_nll] = fmincon(objfun,startvalues_random,[],[],[],[],lb,ub,[],options_fmincon);
    bhat_multi = [bhat_multi;temp_bhat];
    startvalues_multi=[startvalues_multi;startvalues_random];
    nll_multi = [nll_multi;temp_nll]
end
else % if not random
multiplywith=[ones(size(bhat_fmincon,2),size(bhat_fmincon,2))+ diag(ones(1,size(bhat_fmincon,2))), ones(size(bhat_fmincon,2),size(bhat_fmincon,2))+ -1/2 .* diag(ones(1,size(bhat_fmincon,2))) ];    
for i=1:size(multiplywith,2)
    startvalues_random=  startvalues.* multiplywith(:,i)';
    startvalues_random(end)= min(0.9999,max(0.0001,startvalues_random(end)));
    [temp_bhat,temp_nll] = fmincon(objfun,startvalues_random,[],[],[],[],lb,ub,[],options_fmincon);
    bhat_multi = [bhat_multi;temp_bhat];
    startvalues_multi=[startvalues_multi;startvalues_random];
    nll_multi = [nll_multi;temp_nll]
end
end

temp=[nll_multi bhat_multi];
