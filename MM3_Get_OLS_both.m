%% TRAIN THE PCE METAMODELS, TEST THEIR R2
% idea : make a PCE for each cost function under each rain, where the cost
% function is the difference between the simulations and the observation
% coming from the true parameter and the true rain.
clearvars
rng(100,'twister')
cd /Users/icj_ecl/Documents/UQLab_Rel2.1.0/core/
uqlab
output_folder='/Users/icj_ecl/Documents/RESS_review_depo/metamodels/';
data_folder = '/Users/icj_ecl/Documents/RESS_review_depo/'; %/data/raw

%% Set data for the experimental design of the training set
parcel_names = [545,530,480,503,526,527,481,524,525,529,544,539,528,523];
vector_input_indices = [63,68,72,99,71,65]; % Set parameters considered

%% Set data for the experimental design of the training set
R=500; 
R_train = 1:200; % rain realizations for training
R_test = 201:499; % rain realizations used for testing PCE metamodels
N_tot_train = 50; % number of simulations after which the rain changes 

name_train = strcat('LHS_rep5_50_R',int2str(R)); 

%% Load Train set PESHMELBA outputs
cd(data_folder)
lhs_x_train = readNPY(strcat(name_train,'_dim6_sample_final.npy'));
lhs_y_train = readNPY(strcat('Ymoisture_profile_',name_train,'_dim6_YzeronS06_LogN02.npy'));

if (size(lhs_x_train,1) == size(lhs_y_train,1))
    n_lhs_train = size(lhs_x_train,1);
else
    print("error: simulation input and output dimensions incompatible")
end

expdes_x_train = lhs_x_train(1:n_lhs_train, vector_input_indices);
size(expdes_x_train)

%% Load VALIDATION set PESHMELBA outputs
cd(data_folder)
lhs_x_validation = readNPY(strcat('LHS25_100_R500_dim6_sample_final.npy'));
lhs_y_validation = readNPY(strcat('Ymoisture_profile_LHS25_100_R499_dim6_YzeronS06_LogN02.npy'));

if (size(lhs_x_validation,1) == size(lhs_y_validation,1) + 100) % removed last 100 as nominal rain
    n_lhs_validation = size(lhs_y_validation,1);
else
    print("error: simulation input and output dimensions incompatible")
end

expdes_x_validation = lhs_x_validation(1:n_lhs_validation, vector_input_indices);
size(expdes_x_validation)

%% Load true observation set data PESHMELBA, for definition of cost function
cd(data_folder)
% simulations with YzeronS06 rains without rain error
lhs_x_test = readNPY('LHS25_100_R198_dim6_sample_final.npy');
lhs_y_true_rains = readNPY('Ymoisture_profile_LHS25_100_R198_dim6_YzeronS06_hourly.npy');

if (size(lhs_x_test,1) == size(lhs_y_true_rains,1))
    n_lhs_test = size(lhs_x_test,1);
end

% this rain is perturbated with error and used for train simulations.
true_test_idx = 1;
true_rain_idx = 53; % indexed from 0

% true parameter values
x_true  = lhs_x_test(true_test_idx, vector_input_indices);
x_true

% true observation
size(lhs_y_true_rains)
y_true = lhs_y_true_rains(true_rain_idx*100+true_test_idx,:);


%% Reading input parameter MARGINAL distributions, same for all parcels
cd(data_folder)
plot_idx = 4;
load('input_descriptor145.mat')

jj =1;
for ii = vector_input_indices
    InputOpts1.Marginals(jj).Type = deblank(Marginals.Type(ii,:));
    InputOpts1.Marginals(jj).Parameters = Marginals.Parameters{:,ii};
    jj = jj + 1;
end
myInput1 = uq_createInput(InputOpts1);

%% Load OLS basis
beta = 0.005;
str_beta = '005';

cd(output_folder)
A_OLS_union = load(strcat('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.mat'));
cd(output_folder)

%% Calculate the cost function
depths = [0.5 1 2 3 4 5 6 10 15 20 25 30 35 40 45 50 55 65 75 100 150 200 250 300 400];
weighted = 0;
weights = ones(1,length(depths)); % depths

%% Project TRAIN trajectories on the OLS basis
N_train = 50 % number of simulation for PCE fitting

figure
%hold on
for r = 1:200 % number of cost functions to minimize

    % add Jb on the cost function !!!!!!!!!!
    expdes_J_train = ((lhs_y_train(1:n_lhs_train,:) - ones(n_lhs_train,1)*y_true).^2)*weights' + ...
        beta * 1/Marginals.Parameters{:,68}(2).^2 *  (lhs_x_train(1:n_lhs_train,68)- Marginals.Parameters{:,68}(1)).^2+ ...
        beta * 1/Marginals.Parameters{:,71}(2).^2 *  (lhs_x_train(1:n_lhs_train,71)- Marginals.Parameters{:,71}(1)).^2;

    % fit and save PCE MM
    name_file = strcat('pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx),'_Jb',str_beta,'_Ridx', int2str(r));
    cd(data_folder)

    X_temp = expdes_x_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train),:);
    Y_temp = expdes_J_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train));
    
%    calculate_OLS_PCE_parcel(parcel_names(plot_idx),name_file, ...
%        X_temp, Y_temp,14,1, A_OLS_union.a_reduced, ...
%        output_folder);

    for i = 1:6
        subplot(2,3,i)
        hold on
        plot(X_temp(:,i), Y_temp,'bx')
    end

end

%% Project VALIDATION set on the OLS basis
N_valid = 100;      % number of simulation used for PCE fitting
N_tot_valid = 100;  % number of simulations after which the rain changes

for r = 1:499 % number of cost functions

    % add Jb on the cost function !!!!!!!!!!
    expdes_J_validation = ((lhs_y_validation(1:n_lhs_validation,:) - ones(n_lhs_validation,1)*y_true).^2)*weights' + ...
        beta * 1/Marginals.Parameters{:,68}(2).^2 *  (lhs_x_validation(1:n_lhs_validation,68)- Marginals.Parameters{:,68}(1)).^2+ ...
        beta * 1/Marginals.Parameters{:,71}(2).^2 *  (lhs_x_validation(1:n_lhs_validation,71)- Marginals.Parameters{:,71}(1)).^2;
    
    NX_idx = ((r-1)*N_tot_valid+1):((r-1)*N_tot_valid+N_valid)

    % fit and save PCE MM
        name_file = strcat('TEST_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx),'_Jb',str_beta, ...
            '_Ridx', int2str(r));
        cd(data_folder)

    X_temp = expdes_x_validation(NX_idx,:);
    Y_temp = expdes_J_validation(NX_idx);
    
%    calculate_OLS_PCE_parcel(parcel_names(plot_idx),name_file, ...
%        X_temp, Y_temp,14,1, ...
%        A_OLS_union.a_reduced, output_folder);

    for i = 1:6
        subplot(2,3,i)
        hold on
        plot(X_temp(:,i), Y_temp,'ko')
    end

end