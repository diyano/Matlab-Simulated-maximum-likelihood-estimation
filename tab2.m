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


%% Naive
sophisticated=0;

    % Naive(N), With Omega(O), With Sigma(S)
    [~,~,~,~,~,NOS_eu_a,NOS_eu_b] ...
    = solvemodel(bhat1OS(1),[bhat1OS(2);0],bhat1OS(3),bhat1OS(4),phi1OS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), With Omega(O), Without Sigma(s)
    [~,~,~,~,~,NOs_eu_a,NOs_eu_b] ...
    = solvemodel(bhat1Os(1),[bhat1Os(2);0],bhat1Os(3),bhat1Os(4),phi1Os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), Without Omega(o), With Sigma(S)
    [~,~,~,~,~,NoS_eu_a,NoS_eu_b] ...
    = solvemodel(bhat1oS(1),[bhat1oS(2);0],bhat1oS(3),bhat1oS(4),phi1oS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

    % Naive(N), Without Omega(o), Without Sigma(s)
    [~,~,~,~,~,Nos_eu_a,Nos_eu_b] ...
    = solvemodel(bhat1os(1),[bhat1os(2);0],bhat1os(3),bhat1os(4),phi1os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws1,persnr_cont,sophisticated);

%% Sophisticated
sophisticated=1;
 
    % Sophisticated(S), With Omega(O), With Sigma(S)
    [~,~,~,~,~,SOS_eu_a,SOS_eu_b] ...
    = solvemodel(bhat2OS(1),[bhat2OS(2);0],bhat2OS(3),bhat2OS(4),phi2OS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);
 
    % Sophisticated(S), With Omega(O), Without Sigma(s)
    [~,~,~,~,~,SOs_eu_a,SOs_eu_b] ...
    = solvemodel(bhat2Os(1),[bhat2Os(2);0],bhat2Os(3),bhat2Os(4),phi2Os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);
 
    % Sophisticated(S), Without Omega(o), With Sigma(S)
    [~,~,~,~,~,SoS_eu_a,SoS_eu_b] ...
    = solvemodel(bhat2oS(1),[bhat2oS(2);0],bhat2oS(3),bhat2oS(4),phi2oS,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);
 
    % Sophisticated(S), Without Omega(o), Without Sigma(s)
    [~,~,~,~,~,Sos_eu_a,Sos_eu_b] ...
    = solvemodel(bhat2os(1),[bhat2os(2);0],bhat2os(3),bhat2os(4),phi2os,a1,a2,b1,b2,typeq,...
    X_beta,X_delta,X_sigma,X_omega,standardnormaldraws2,persnr_cont,sophisticated);

%% find the averages for each question

%typeq_extended= 8127*100x1
typeq_extended = reshape( repmat(typeq,1,numberdraws),[],1);

%NOS_eu_a_extended= 8127*100x1
NOS_eu_a_extended= reshape(NOS_eu_a,[],1);
NOS_eu_b_extended= reshape(NOS_eu_b,[],1);
NOs_eu_a_extended= reshape(NOs_eu_a,[],1);
NOs_eu_b_extended= reshape(NOs_eu_b,[],1);
NoS_eu_a_extended= reshape(NoS_eu_a,[],1);
NoS_eu_b_extended= reshape(NoS_eu_b,[],1);
Nos_eu_a_extended= reshape(Nos_eu_a,[],1);
Nos_eu_b_extended= reshape(Nos_eu_b,[],1);
SOS_eu_a_extended= reshape(SOS_eu_a,[],1);
SOS_eu_b_extended= reshape(SOS_eu_b,[],1);
SOs_eu_a_extended= reshape(SOs_eu_a,[],1);
SOs_eu_b_extended= reshape(SOs_eu_b,[],1);
SoS_eu_a_extended= reshape(SoS_eu_a,[],1);
SoS_eu_b_extended= reshape(SoS_eu_b,[],1);
Sos_eu_a_extended= reshape(Sos_eu_a,[],1);
Sos_eu_b_extended= reshape(Sos_eu_b,[],1);

%temp_NOS_eu_a_extended= 2709X100
temp_NOS_eu_a_extended=reshape(NOS_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOS_eu_b_extended=reshape(NOS_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOs_eu_a_extended=reshape(NOs_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NOs_eu_b_extended=reshape(NOs_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoS_eu_a_extended=reshape(NoS_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_NoS_eu_b_extended=reshape(NoS_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_Nos_eu_a_extended=reshape(Nos_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_Nos_eu_b_extended=reshape(Nos_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOS_eu_a_extended=reshape(SOS_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOS_eu_b_extended=reshape(SOS_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOs_eu_a_extended=reshape(SOs_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SOs_eu_b_extended=reshape(SOs_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoS_eu_a_extended=reshape(SoS_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_SoS_eu_b_extended=reshape(SoS_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_Sos_eu_a_extended=reshape(Sos_eu_a_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);
temp_Sos_eu_b_extended=reshape(Sos_eu_b_extended(typeq_extended==4|typeq_extended==5|typeq_extended==6),[],numberdraws);

%plot_NOS_eu_a = 21 (means taken over draws and then subjects)
plot_NOS_eu_a = mean(mean(reshape(temp_NOS_eu_a_extended,21,[],numberdraws),3),2);
plot_NOS_eu_b = mean(mean(reshape(temp_NOS_eu_b_extended,21,[],numberdraws),3),2);
plot_NOs_eu_a = mean(mean(reshape(temp_NOs_eu_a_extended,21,[],numberdraws),3),2);
plot_NOs_eu_b = mean(mean(reshape(temp_NOs_eu_b_extended,21,[],numberdraws),3),2);
plot_NoS_eu_a = mean(mean(reshape(temp_NoS_eu_a_extended,21,[],numberdraws),3),2);
plot_NoS_eu_b = mean(mean(reshape(temp_NoS_eu_b_extended,21,[],numberdraws),3),2);
plot_Nos_eu_a = mean(mean(reshape(temp_Nos_eu_a_extended,21,[],numberdraws),3),2);
plot_Nos_eu_b = mean(mean(reshape(temp_Nos_eu_b_extended,21,[],numberdraws),3),2);

plot_SOS_eu_a = mean(mean(reshape(temp_SOS_eu_a_extended,21,[],numberdraws),3),2);
plot_SOS_eu_b = mean(mean(reshape(temp_SOS_eu_b_extended,21,[],numberdraws),3),2);
plot_SOs_eu_a = mean(mean(reshape(temp_SOs_eu_a_extended,21,[],numberdraws),3),2);
plot_SOs_eu_b = mean(mean(reshape(temp_SOs_eu_b_extended,21,[],numberdraws),3),2);
plot_SoS_eu_a = mean(mean(reshape(temp_SoS_eu_a_extended,21,[],numberdraws),3),2);
plot_SoS_eu_b = mean(mean(reshape(temp_SoS_eu_b_extended,21,[],numberdraws),3),2);
plot_Sos_eu_a = mean(mean(reshape(temp_Sos_eu_a_extended,21,[],numberdraws),3),2);
plot_Sos_eu_b = mean(mean(reshape(temp_Sos_eu_b_extended,21,[],numberdraws),3),2);

%% discount by beta delta
delta1=mean(max(0,1-par_draws1(:,2).^2),1);
beta1=mean(par_draws1(:,1),1);
delta2=mean(max(0,1-par_draws2(:,2).^2),1);
beta2=mean(par_draws2(:,1),1);
plot_NOS_eu_a = plot_NOS_eu_a * delta1 * beta1; 
plot_NOS_eu_b = plot_NOS_eu_b * delta1 * beta1; 
plot_NOs_eu_a = plot_NOs_eu_a * delta1 * beta1; 
plot_NOs_eu_b = plot_NOs_eu_b * delta1 * beta1; 
plot_NoS_eu_a = plot_NoS_eu_a * delta1 * beta1; 
plot_NoS_eu_b = plot_NoS_eu_b * delta1 * beta1; 
plot_Nos_eu_a = plot_Nos_eu_a * delta1 * beta1; 
plot_Nos_eu_b = plot_Nos_eu_b * delta1 * beta1; 

plot_SOS_eu_a = plot_SOS_eu_a * delta2 * beta2; 
plot_SOS_eu_b = plot_SOS_eu_b * delta2 * beta2; 
plot_SOs_eu_a = plot_SOs_eu_a * delta2 * beta2; 
plot_SOs_eu_b = plot_SOs_eu_b * delta2 * beta2; 
plot_SoS_eu_a = plot_SoS_eu_a * delta2 * beta2; 
plot_SoS_eu_b = plot_SoS_eu_b * delta2 * beta2;
plot_Sos_eu_a = plot_Sos_eu_a * delta2 * beta2;
plot_Sos_eu_b = plot_Sos_eu_b * delta2 * beta2;

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

input.data = [payments plot_NOS_eu_a plot_NOs_eu_a plot_NOS_eu_b plot_NOs_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A','$\sigma = 0$, A','$\sigma \neq 0$, B','$\sigma = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\omega \neq 0$';
% LaTex table label:
input.tableLabel = 'naivewithomega';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%No
fprintf('%% Table 2, Naive without omega \n')
input.data = [payments plot_NoS_eu_a plot_Nos_eu_a plot_NoS_eu_b plot_Nos_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A','$\sigma = 0$, A','$\sigma \neq 0$, B','$\sigma = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\omega = 0$';
% LaTex table label:
input.tableLabel = 'naivewithoutomega';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%SO
fprintf('%% Table 2, Sophisticated with omega \n')
input.data = [payments plot_SOS_eu_a plot_SOs_eu_a plot_SOS_eu_b plot_SOs_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A','$\sigma = 0$, A','$\sigma \neq 0$, B','$\sigma = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\omega \neq 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithomega';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%So
fprintf('%% Table 2, Sophisticated without omega \n')
input.data = [payments plot_SoS_eu_a plot_Sos_eu_a plot_SoS_eu_b plot_Sos_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\sigma \neq 0$, A','$\sigma = 0$, A','$\sigma \neq 0$, B','$\sigma = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\omega = 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithoutomega';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

%% tables for 4 categories NS , Ns, SS, Ss

% NS

fprintf('%% Table 2, Naive with sigma \n')

input.data = [payments plot_NOS_eu_a plot_NoS_eu_a plot_NOS_eu_b plot_NoS_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A','$\omega = 0$, A','$\omega \neq 0$, B','$\omega = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\sigma \neq 0$';
% LaTex table label:
input.tableLabel = 'naivewithsigma';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% Ns

fprintf('%% Table 2, Naive without sigma \n')

input.data = [payments plot_NOs_eu_a plot_Nos_eu_a plot_NOs_eu_b plot_Nos_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A','$\omega = 0$, A','$\omega \neq 0$, B','$\omega = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Naive, $\sigma = 0$';
% LaTex table label:
input.tableLabel = 'naivewithoutsigma';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% SS

fprintf('%% Table 2, Sophisticated with sigma \n')

input.data = [payments plot_SOS_eu_a plot_SoS_eu_a plot_SOS_eu_b plot_SoS_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A','$\omega = 0$, A','$\omega \neq 0$, B','$\omega = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\sigma \neq 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithsigma';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

% Ss

fprintf('%% Table 2, Sophisticated without sigma \n')

input.data = [payments plot_SOs_eu_a plot_Sos_eu_a plot_SOs_eu_b plot_Sos_eu_b];
input.tableColLabels = { 'Type'  'Question' '$f_2$' '$f_4$' '$g_2$' '$g_4$' '$\omega \neq 0$, A','$\omega = 0$, A','$\omega \neq 0$, B','$\omega = 0$, B'};
input.tableRowLabels = { ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
input.transposeTable = 0;
% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat = {'%.0f',2,'%.2f',8}; % two digits precision for first two columns, two digit for the last two
% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';
% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sophisticated, $\sigma = 0$';
% LaTex table label:
input.tableLabel = 'sophisticatedwithoutsigma';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% call latexTable:
latex = latexTable(input);

fprintf('%% END OF COPY AND PASTE \n')

%% figures for 4 categories NO , No, SO, So

%NO
close all
subplot(1,2,1) 
bar(1:21,[plot_NOS_eu_a plot_NOs_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOS_eu_a plot_NOs_eu_a plot_NOS_eu_b plot_NOs_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOS_eu_a plot_NOs_eu_a plot_NOS_eu_b plot_NOs_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_NOS_eu_b plot_NOs_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOS_eu_a plot_NOs_eu_a plot_NOS_eu_b plot_NOs_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOS_eu_a plot_NOs_eu_a plot_NOS_eu_b plot_NOs_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Naive \omega \neq 0')  

print('-dpdf','../output/matlab_table2_NbigO.pdf')

%No
close all
subplot(1,2,1) 
bar(1:21,[plot_NoS_eu_a plot_Nos_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NoS_eu_a plot_Nos_eu_a plot_NoS_eu_b plot_Nos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NoS_eu_a plot_Nos_eu_a plot_NoS_eu_b plot_Nos_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \omega = 0')  

subplot(1,2,2) 
bar(1:21,[plot_NoS_eu_b plot_Nos_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NoS_eu_a plot_Nos_eu_a plot_NoS_eu_b plot_Nos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NoS_eu_a plot_Nos_eu_a plot_NoS_eu_b plot_Nos_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Naive \omega = 0')  

print('-dpdf','../output/matlab_table2_Nsmallo.pdf')

%SO
close all
subplot(1,2,1) 
bar(1:21,[plot_SOS_eu_a plot_SOs_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOS_eu_a plot_SOs_eu_a plot_SOS_eu_b plot_SOs_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOS_eu_a plot_SOs_eu_a plot_SOS_eu_b plot_SOs_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_SOS_eu_b plot_SOs_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOS_eu_a plot_SOs_eu_a plot_SOS_eu_b plot_SOs_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOS_eu_a plot_SOs_eu_a plot_SOS_eu_b plot_SOs_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Sophisticated \omega \neq 0')  

print('-dpdf','../output/matlab_table2_SbigO.pdf')

%So
close all
subplot(1,2,1) 
bar(1:21,[plot_SoS_eu_a plot_Sos_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SoS_eu_a plot_Sos_eu_a plot_SoS_eu_b plot_Sos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SoS_eu_a plot_Sos_eu_a plot_SoS_eu_b plot_Sos_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \omega = 0')  

subplot(1,2,2) 
bar(1:21,[plot_SoS_eu_b plot_Sos_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SoS_eu_a plot_Sos_eu_a plot_SoS_eu_b plot_Sos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SoS_eu_a plot_Sos_eu_a plot_SoS_eu_b plot_Sos_eu_b]))]);
hleg1 = legend('\sigma \neq 0','\sigma = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Sophisticated \omega = 0')  

print('-dpdf','../output/matlab_table2_Ssmallo.pdf')

%% figures for 4 categories NS , Ns, SS, Ss
% NS
close all
subplot(1,2,1) 
bar(1:21,[plot_NOS_eu_a plot_NoS_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOS_eu_a plot_NoS_eu_a plot_NOS_eu_b plot_NoS_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOS_eu_a plot_NoS_eu_a plot_NOS_eu_b plot_NoS_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_NOS_eu_b plot_NoS_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOS_eu_a plot_NoS_eu_a plot_NOS_eu_b plot_NoS_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOS_eu_a plot_NoS_eu_a plot_NOS_eu_b plot_NoS_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Naive \sigma \neq 0')  

print('-dpdf','../output/matlab_table2_NbigS.pdf')

% Ns
close all
subplot(1,2,1) 
bar(1:21,[plot_NOs_eu_a plot_Nos_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOs_eu_a plot_Nos_eu_a plot_NOs_eu_b plot_Nos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOs_eu_a plot_Nos_eu_a plot_NOs_eu_b plot_Nos_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Naive \sigma = 0')  

subplot(1,2,2) 
bar(1:21,[plot_NOs_eu_b plot_Nos_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_NOs_eu_a plot_Nos_eu_a plot_NOs_eu_b plot_Nos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_NOs_eu_a plot_Nos_eu_a plot_NOs_eu_b plot_Nos_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Naive \sigma = 0')  

print('-dpdf','../output/matlab_table2_NsmallS.pdf')

% SS
close all
subplot(1,2,1) 
bar(1:21,[plot_SOS_eu_a plot_SoS_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOS_eu_a plot_SoS_eu_a plot_SOS_eu_b plot_SoS_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOS_eu_a plot_SoS_eu_a plot_SOS_eu_b plot_SoS_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma \neq 0')  

subplot(1,2,2) 
bar(1:21,[plot_SOS_eu_b plot_SoS_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOS_eu_a plot_SoS_eu_a plot_SOS_eu_b plot_SoS_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOS_eu_a plot_SoS_eu_a plot_SOS_eu_b plot_SoS_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Sophisticated \sigma \neq 0')  

print('-dpdf','../output/matlab_table2_SbigS.pdf')

% Ss
close all
subplot(1,2,1) 
bar(1:21,[plot_SOs_eu_a plot_Sos_eu_a]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOs_eu_a plot_Sos_eu_a plot_SOs_eu_b plot_Sos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOs_eu_a plot_Sos_eu_a plot_SOs_eu_b plot_Sos_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option A, Sophisticated \sigma = 0')  

subplot(1,2,2) 
bar(1:21,[plot_SOs_eu_b plot_Sos_eu_b]);
barmap=[207/256 181/256 59/256; 30/256 144/256 255/256]; 
colormap(barmap);
axis([0 22 0 1.2*max(max([plot_SOs_eu_a plot_Sos_eu_a plot_SOs_eu_b plot_Sos_eu_b]))]);
set(gca,'XTick',1:21)
set(gca,'YTick',[0:2:1.2*max(max([plot_SOs_eu_a plot_Sos_eu_a plot_SOs_eu_b plot_Sos_eu_b]))]);
hleg1 = legend('\omega \neq 0','\omega = 0');
xlabel('Question','FontSize',fontsize);
ylabel('Probability','FontSize',fontsize);
set(hleg1,'Location','NorthEast');
title('Option B, Sophisticated \sigma = 0')  

print('-dpdf','../output/matlab_table2_Ssmalls.pdf')


%% value of flexibility naive
mean((plot_NOS_eu_a - plot_NOs_eu_a) + (plot_NOS_eu_a - plot_NoS_eu_a))
mean((plot_NOS_eu_b - plot_NOs_eu_b) + (plot_NOS_eu_b - plot_NoS_eu_b))
mean((plot_NOS_eu_a - plot_NOs_eu_a) + (plot_NOS_eu_a - plot_NoS_eu_a)) / mean(plot_NOS_eu_a)

% value of flexibility sophisticated
mean((plot_SOS_eu_a - plot_SOs_eu_a) + (plot_SOS_eu_a - plot_SoS_eu_a))
mean((plot_SOS_eu_b - plot_SOs_eu_b) + (plot_SOS_eu_b - plot_SoS_eu_b))
mean((plot_SOS_eu_a - plot_SOs_eu_a) + (plot_SOS_eu_a - plot_SoS_eu_a)) / mean(plot_SOS_eu_a)