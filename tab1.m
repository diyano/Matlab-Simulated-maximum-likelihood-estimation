close all
%arrange the results
results=[ bhat1 se_bhat1 ; phi1 se_phi1 ; bhat2 se_bhat2 ; phi2 se_phi2; bhat(size(bhat,1)) se(size(se,1))] ;

results(results==0)=NaN;

fprintf('%% \n %% \n %% \n %% \n %% \n')
fprintf('%% COPY AND PASTE STARTING FROM HERE \n') % this marks the part which you should copy from the console
fprintf('%% Table 1, Estimation results \n')
input.data = results;
input.tableColLabels = { 'Estimates'  'Std. error' };
input.tableRowLabels = { '$\mu_{\beta^{*}_{N}}$ ' '$\mu_{\delta^{*}_{N}}$ ' '$\mu_{\sigma^{*}_{N}}$' '$\mu_{\omega^{*}_{N}}$ '...
'$\rho_{\beta^{*}_{N},\beta^{*}_{N}}$' '$\rho_{\beta^{*}_{N},\delta^{*}_{N}}$' '$\rho_{\delta^{*}_{N},\delta^{*}_{N}}$' '$\rho_{\sigma^{*}_{N},\beta^{*}_{N}}$' '$\rho_{\sigma^{*}_{N},\delta^{*}_{N}}$' '$\rho_{\sigma^{*}_{N},\sigma^{*}_{N}}$' '$\rho_{\omega^{*}_{N},\beta^{*}_{N}}$' '$\rho_{\omega^{*}_{N},\delta^{*}_{N}}$' '$\rho_{\omega^{*}_{N},\sigma^{*}_{N}}$' '$\rho_{\omega^{*}_{N},\omega^{*}_{N}}$'...   
'$\mu_{\beta^{*}_{S}}$ ' '$\mu_{\delta^{*}_{S}}$ ' '$\mu_{\sigma^{*}_{S}}$' '$\mu_{\omega^{*}_{S}}$ '...
'$\rho_{\beta^{*}_{S},\beta^{*}_{S}}$' '$\rho_{\beta^{*}_{S},\delta^{*}_{S}}$' '$\rho_{\delta^{*}_{S},\delta^{*}_{S}}$' '$\rho_{\sigma^{*}_{S},\beta^{*}_{S}}$' '$\rho_{\sigma^{*}_{S},\delta^{*}_{S}}$' '$\rho_{\sigma^{*}_{S},\sigma^{*}_{S}}$' '$\rho_{\omega^{*}_{S},\beta^{*}_{S}}$' '$\rho_{\omega^{*}_{S},\delta^{*}_{S}}$' '$\rho_{\omega^{*}_{S},\sigma^{*}_{S}}$' '$\rho_{\omega^{*}_{S},\omega^{*}_{S}}$'...
'$\pi$'};
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
input.tableCaption = 'Estimation results';
% LaTex table label:
input.tableLabel = 'estresults';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);
fprintf('%% PLEASE WAIT \n')
