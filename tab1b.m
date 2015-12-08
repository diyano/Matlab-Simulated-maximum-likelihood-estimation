%% plot distribution random coefficients
close all
% naive
draws_beta1 = par_draws1(:,1);
draws_delta1 = max(0,1-par_draws1(:,2).^2);
draws_sigma1 = max(1e-1,min(par_draws1(:,3).^2,10));
draws_omega1 = par_draws1(:,4);

% sophisticated
draws_beta2 = par_draws2(:,1);
draws_delta2 = max(0,1-par_draws2(:,2).^2);
draws_sigma2 = max(1e-1,min(par_draws2(:,3).^2,10));
draws_omega2 = par_draws2(:,4);

paren = @(x, varargin) x(varargin{:});

original=[bhat1;phi1;bhat2;phi2;bhat_fmincon(end)];

transformed= [ mean(draws_beta1); mean(draws_delta1); mean(draws_sigma1); mean(draws_omega1); ...
std(draws_beta1); paren(corrcoef(draws_beta1,draws_delta1),1,2); std(draws_delta1);  paren(corrcoef(draws_sigma1,draws_beta1),1,2);  paren(corrcoef(draws_sigma1,draws_delta1),1,2); std(draws_sigma1); paren(corrcoef(draws_omega1,draws_beta1),1,2); paren(corrcoef(draws_omega1,draws_delta1),1,2); ; paren(corrcoef(draws_omega1,draws_sigma1),1,2); std(draws_omega1); ...   
mean(draws_beta2); mean(draws_delta2); mean(draws_sigma2); mean(draws_omega2); ... 
std(draws_beta2); paren(corrcoef(draws_beta2,draws_delta2),1,2); std(draws_delta2);  paren(corrcoef(draws_sigma2,draws_beta2),1,2);  paren(corrcoef(draws_sigma2,draws_delta2),1,2); std(draws_sigma2); paren(corrcoef(draws_omega2,draws_beta2),1,2); paren(corrcoef(draws_omega2,draws_delta2),1,2); ; paren(corrcoef(draws_omega2,draws_sigma2),1,2); std(draws_omega2) ...
];
transformed(transformed==0)=NaN;
transformed(original==0)=NaN;
transformed=[NaN(28,1) , transformed];
%transformed=[5, 6 ; 8,9];
% paren(corrcoef(draws_beta1,draws_delta1),1,2)

input.data = transformed;
input.tableColLabels = { '' '' };
% input.tableRowLabels = {'label1' 'label2' };
input.tableRowLabels = { '$\mu_{\beta_{N}}$ '  '$\mu_{\delta_{N}}$ ' '$\mu_{\sigma_{N}}$' '$\mu_{\omega_{N}}$ '...
'$\rho_{\beta_{N},\beta_{N}}$' '$\rho_{\beta_{N},\delta_{N}}$' '$\rho_{\delta_{N},\delta_{N}}$' '$\rho_{\sigma_{N},\beta_{N}}$' '$\rho_{\sigma_{N},\delta_{N}}$' '$\rho_{\sigma_{N},\sigma_{N}}$' '$\rho_{\omega_{N},\beta_{N}}$' '$\rho_{\omega_{N},\delta_{N}}$' '$\rho_{\omega_{N},\sigma_{N}}$' '$\rho_{\omega_{N},\omega_{N}}$'...   
'$\mu_{\beta_{S}}$ ' '$\mu_{\delta_{S}}$ ' '$\mu_{\sigma_{S}}$' '$\mu_{\omega_{S}}$ '...
'$\rho_{\beta_{S},\beta_{S}}$' '$\rho_{\beta_{S},\delta_{S}}$' '$\rho_{\delta_{S},\delta_{S}}$' '$\rho_{\sigma_{S},\beta_{S}}$' '$\rho_{\sigma_{S},\delta_{S}}$' '$\rho_{\sigma_{S},\sigma_{S}}$' '$\rho_{\omega_{S},\beta_{S}}$' '$\rho_{\omega_{S},\delta_{S}}$' '$\rho_{\omega_{S},\sigma_{S}}$' '$\rho_{\omega_{S},\omega_{S}}$'};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.3f',2}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Transformed parameters';
% LaTex table label:
input.tableLabel = 'transformed';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

