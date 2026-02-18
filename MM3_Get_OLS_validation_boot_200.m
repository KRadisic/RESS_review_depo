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

%% Additional simulations read

cd(data_folder)
lhs_x_add = readNPY(strcat('LHS110_200_R300_dim6_sample_final.npy'));
lhs_y_add = readNPY(strcat('Ymoisture_profile_LHS91_200_R200_499_dim6_YzeronS06_LogN02.npy'));

if (size(lhs_x_add,1) == size(lhs_y_add,1))
    n_lhs_add = size(lhs_y_add,1);
else
    print("error: simulation input and output dimensions incompatible")
end

%% Additional simulations concatenate

lhs_y_validation_add = [];
expdes_x_validation_add = [];
lhs_x_validation_add = [];

% r_idx = 20
% nidx = 95
% lhs_y_validation(nidx + (100*(200+r_idx-1-1)),1:25) %corresponds to sample(110*(317-200)+97-90-1) new simus
% lhs_y_add(nidx - 90 + (r_idx-1)*110, 1:25)
% 

for r_idx = 1:299 % check the last rain ??
    idxs_old = (1 + 100*(200+r_idx-1-1)) : (100*(200+r_idx-1));
    idxs_new = (11 + (r_idx-1)*110) : (110 + (r_idx-1)*110);

    lhs_y_validation_add = [lhs_y_validation_add; lhs_y_validation(idxs_old,1:25)];
    lhs_y_validation_add = [lhs_y_validation_add; lhs_y_add(idxs_new,1:25)];

    expdes_x_validation_add = [expdes_x_validation_add; expdes_x_validation(idxs_old,1:6)];
    expdes_x_validation_add = [expdes_x_validation_add; lhs_x_add(idxs_new,vector_input_indices)];

    lhs_x_validation_add = [lhs_x_validation_add; lhs_x_validation(idxs_old,:)];
    lhs_x_validation_add = [lhs_x_validation_add; lhs_x_add(idxs_new,:)];

end

size(expdes_x_validation_add)
size(lhs_y_validation_add)
size(lhs_x_validation_add)

%% Project VALIDATION set on the OLS basis
%Nx_valid = 100;      % number of simulation used for PCE fitting
N_valid_boot = 8;  % number of bootstraps (how many subsets of 50 points out of 100 we take) 
N_tot_valid = 200;  % number of simulations after which the rain changes

n_lhs_validation_add = 59800;

for Nx_valid = [50, 100]     % number of simulation used for PCE fitting
    figure
    for i_boot = 1:N_valid_boot % idx of bootstrap
        NX_idx_boot = randperm(N_tot_valid, Nx_valid) % randomly choose Nx_valid indices without repetition from N_tot_valid
        % indexing starts from 1, goes to 100
        
        for r = 1:299 % number of cost functions
            % add Jb on the cost function !!!!!!!!!!
            expdes_J_validation = ((lhs_y_validation_add(1:n_lhs_validation_add,:) - ones(n_lhs_validation_add,1)*y_true).^2)*weights' + ...
                beta * 1/Marginals.Parameters{:,68}(2).^2 *  (lhs_x_validation_add(1:n_lhs_validation_add,68)- Marginals.Parameters{:,68}(1)).^2+ ...
                beta * 1/Marginals.Parameters{:,71}(2).^2 *  (lhs_x_validation_add(1:n_lhs_validation_add,71)- Marginals.Parameters{:,71}(1)).^2;
            
            NX_idx = ((r-1) * N_tot_valid + NX_idx_boot)
        
            % fit and save PCE MM
            name_file = strcat('TEST_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx),'_Jb',str_beta, ...
                '_Ridx', int2str(r+200), '_boot', int2str(i_boot), '_of', int2str(N_valid_boot), 'Nx_valid', int2str(Nx_valid));
        
            X_temp = expdes_x_validation_add(NX_idx,:);
            Y_temp = expdes_J_validation(NX_idx);
            
            cd(data_folder)
            calculate_OLS_PCE_parcel(parcel_names(plot_idx),name_file, ...
            X_temp, Y_temp,14,1, A_OLS_union.a_reduced, output_folder);
            
            for i = 1:6
                subplot(2,3,i)
                hold on
                plot(X_temp(:,i), Y_temp,'bx')
            end
        end
    end
end