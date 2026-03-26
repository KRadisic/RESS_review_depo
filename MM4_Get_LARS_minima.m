%% Minimize train trjectories
clearvars
rng(100,'twister')
cd '/home/kradisic/Documents/my_UQLab_folder/UQLab_Rel2.2.0/core' 

uqlab
output_folder='/home/kradisic/Documents/RESS_review_depo/metamodels';
data_folder = '/home/kradisic/Documents/RESS_review_depo/';

%% Set data
R=500; 
name_train = strcat('LHS_rep5_50_R',int2str(R)); 

%% Set parameters considered
vector_input_indices = [63,68,72,99,71,65];
N_train = 50; % number of simulations after which the rain changes 

N_tot_train =50;
weighted = 0;
true_rain_idx = 53;

%% Minimisation with the PCE LARS, do we find theta star ?
X_min_LARS = zeros(R,length(vector_input_indices));
f_min_LARS = zeros(R,1);

beta = 0.005;
str_beta = '005';

for rain_idx = 1:200 % loop on train rains
    cd(output_folder)
    load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx),'_Jb',str_beta,'_Ridx', int2str(rain_idx), 'N_train', int2str(N_train),'_503'));
    evaluate_pce_handle = @(x) uq_evalModel(myPCE_LARS,x);

    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [0.229577, 0.1237539, 0.2183556, 0.2332047, 0.01391467, -3.6]; % fix quantiles
    ub = [0.434923, 0.2344461, 0.4136644, 0.4417953, 0.10844333, -2]; % fix quantiles
    x0 = [0.31,0.20,0.34,0.28,0.07,-3.6];
    [xmin,fval,exitflag,output] = fmincon(evaluate_pce_handle,x0,A,b,Aeq,beq,lb,ub);  
    X_min_LARS(rain_idx,:) = xmin;
    f_min_LARS(rain_idx) = fval;
end
% quantiles 0.999 and 0.001
%   print(c(qnorm(0.001, mean = prior_distrib_param$mean, sd = prior_distrib_param$sd),
%          qnorm(0.999, mean = prior_distrib_param$mean, sd = prior_distrib_param$sd)))
% [1] 0.229577 0.434923
% [1] 0.1237539 0.2344461
% [1] 0.2183556 0.4136644
% [1] 0.2332047 0.4417953
% [1] 0.01391467 0.10844333

%% quick plot histogram
figure
vector_input_names = ["thetas9","mn10","thetas10","thetas13","thetar10","hg10"];

for i = 1:6
subplot(2,3,i)
%hist(Sobol_Total_LARS(i,:));
hist(X_min_LARS(:,i))
title(vector_input_names{i})
xlim([lb(i), ub(i)])
end

%% save conditional minimizers and the value at the minimum
cd(output_folder)
X_min_LARS_heading = [[vector_input_names,"fcondmin"]; [X_min_LARS,f_min_LARS]];

writematrix(X_min_LARS_heading,strcat('conditional_minimizers_LARS_errorLogN02_Jb_beta',str_beta,name_train,'N_train', int2str(N_train),'_hourly_BDDquantiles.csv'));
