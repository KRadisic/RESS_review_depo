%% TRAIN THE PCE METAMODELS, TEST THEIR R2
% idea : make a PCE for each cost function under each rain, where the cost
% function is the difference between the simulations and the observation
% coming from the true parameter and the true rain.

% read train set peshmelba simulations
% fit PCE LARS MM on the first 200

clearvars
rng(100,'twister')
cd '/home/katarina.radisic/UQLab_Rel2.0.0/core/'
uqlab
output_folder='/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels';
code_folder = '/home/katarina.radisic/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/2_Robustness_scripts_Abasis';
data_folder = '/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw';

%% Set data for the experimental design of the training set
parcel_names = [545,530,480,503,526,527,481,524,525,529,544,539,528,523];
R=500; 
R_train = 1:200; % rain realizations for training
R_test = 201:499; % rain realizations used tof testing PCE metamodels
N_tot_train = 50; % number of simulations after which the rain changes 

name_train = strcat('LHS_rep5_50_R',int2str(R)); 
vector_input_indices = [63,68,72,99,71,65]; % Set parameters considered

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

%% Calculate the cost function
depths = [0.5 1 2 3 4 5 6 10 15 20 25 30 35 40 45 50 55 65 75 100 150 200 250 300 400];
weighted = 0;
weights = ones(1,length(depths)); % depths

%% Reading input parameter marginal distributions, same for all parcels
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

%% Train all the PCE metamodels

beta = 0.005;
str_beta = '005';
%%
for N_train = 30:5:50 % number of simulation for PCE fitting
    for r = 1:200 %1:198 % number of cost functions to minimize

        % add Jb on the cost function !!!!!!!!!!
        expdes_J_train = ((lhs_y_train(1:n_lhs_train,:) - ones(n_lhs_train,1)*y_true).^2)*weights' + ...
            beta * 1/Marginals.Parameters{:,68}(2).^2 *  (lhs_x_train(1:n_lhs_train,68)- Marginals.Parameters{:,68}(1)).^2+ ...
            beta * 1/Marginals.Parameters{:,71}(2).^2 *  (lhs_x_train(1:n_lhs_train,71)- Marginals.Parameters{:,71}(1)).^2;
        % fit and save PCE MM
            name_file = strcat('pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx), ...
                '_Jb',str_beta,'_Ridx', int2str(r), 'N_train', int2str(N_train));
            cd(code_folder)
            calculate_LARS_PCE_parcel(parcel_names(plot_idx),name_file, ...
                expdes_x_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train),:), ...
                expdes_J_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train)),14,1, ...
                output_folder);
    end
end

%% Get R2 on train set of all the PCE metamodels
coeff_R2_rains = []
figure

for r = 1:200%1:198 % number of cost functions to minimize
    cd(output_folder)
    load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx), ...
        '_Jb',str_beta,'_Ridx', int2str(r),'_503.mat'))
    % get R2
    JLARS = uq_evalModel(myPCE_LARS,myPCE_LARS.ExpDesign.X);
    Jval = myPCE_LARS.ExpDesign.Y;

    plot(Jval,JLARS,'bo')
    hold on
    SSE = sum((Jval - JLARS).^2); 
    SST = sum((Jval-mean(Jval)).^2);
    coeff_R2 = 1 - SSE/SST
    coeff_R2_rains = [coeff_R2_rains,coeff_R2];
end
%%
find(coeff_R2_rains<0.95)
coeff_R2_rains(coeff_R2_rains<0.95)


%% Get R2 on metamodels (validate on the test set available)
all_Q2_rains = []; N_tot_test = 100;
N_valid = 100;      % number of simulation used for PCE fitting
N_tot_valid = 100;  % number of simulations after which the rain changes

figure
for N_train = 30:5:50
    coeff_Q2_rains = [];
    R_temp = 200; 
    
    for rain_idx = 1:R_temp  
        % adapt cost function with Jb
         expdes_J_validation = ((lhs_y_validation(1:n_lhs_validation,:) - ones(n_lhs_validation,1)*y_true).^2)*weights' + ...
            beta * 1/Marginals.Parameters{:,68}(2).^2 *  (lhs_x_validation(1:n_lhs_validation,68)- Marginals.Parameters{:,68}(1)).^2+ ...
            beta * 1/Marginals.Parameters{:,71}(2).^2 *  (lhs_x_validation(1:n_lhs_validation,71)- Marginals.Parameters{:,71}(1)).^2;
        
        cd(output_folder)
        load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx), ...
            '_Jb',str_beta,'_Ridx', int2str(rain_idx),'N_train', int2str(N_train),'_503.mat'))
        
        Xval = expdes_x_validation((2+(rain_idx-1)*N_tot_test):(rain_idx*N_tot_test),:); 
        Jval = expdes_J_validation((2+(rain_idx-1)*N_tot_test):(rain_idx*N_tot_test));
        JLARS = uq_evalModel(myPCE_LARS,Xval);
        %figure
        hold on
        plot(Jval,JLARS,'bo')
        SSE = sum((Jval - JLARS).^2); 
        SST = sum((Jval-mean(Jval)).^2);
        coeff_Q2 = 1 - SSE/SST
        coeff_Q2_rains = [coeff_Q2_rains,coeff_Q2];

        % save test vs metamodel for each rain
        cd(output_folder)
        writematrix([Jval,JLARS],strcat('train_vs_test_PCE_rain',int2str(rain_idx),'.csv'))

    end
    all_Q2_rains = [all_Q2_rains, coeff_Q2_rains'];
end
plot(xlim,ylim,'-k', 'LineWidth',1)

%% Save Q2 for pretty plotting
cd(output_folder)
writematrix(all_Q2_rains,strcat('Q2_LARS_pesh_profmoist_Jpce_errorLogN02_truerain',int2str(true_rain_idx), ...
            '_Jb',str_beta,'_Ntrain_rep5_50_dim6_moistureprofile_YzeronS06_hourly.csv'))
