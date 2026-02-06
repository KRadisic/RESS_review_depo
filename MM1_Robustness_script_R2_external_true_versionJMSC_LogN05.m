%% TRAIN THE PCE METAMODELS, TEST THEIR R2
% idea : make a PCE for each cost function under each rain, where the cost
% function is the difference between the simulations and the observation
% coming from the true parameter and the true rain.
clearvars
rng(100,'twister')
cd '/home/katarina.radisic/UQLab_Rel2.0.0/core/'
uqlab
output_folder='/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels';
code_folder = '/home/katarina.radisic/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/2_Robustness_scripts_Abasis';
data_folder = '/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw';

%% Set data
R=500; 
name_train = strcat('LHS_rep5_50_R',int2str(R)); 

%% Set parameters considered
vector_input_indices = [63,68,72,99,71,65];
%straring from 0 the indices of zA and zB;  0 - 0, 1 - 53, 2 - 54
% mn = 68, thetar = 71

%% Fit PCE LARS
N_tot_test = 100; % number of simulations after which the rain changes 
% fix the index at which we consider the "true" rain observation is coming
% from (taken from the "test" set)
% one simulation of the test set is used as the true yobs for
% calculating the cost function. (and for nothing else in this case)

%% Load Train set PESHMELBA
cd(data_folder)
lhs_x_train = readNPY(strcat(name_train,'_dim6_sample_final.npy'));
lhs_y_train = readNPY(strcat('Ymoisture_profile_',name_train,'_dim6_YzeronS06_LogN.npy'));

if (size(lhs_x_train,1) == size(lhs_y_train,1))
    n_lhs_train = size(lhs_x_train,1);
else
    print("error: simulation input and output dimensions incompatible")
end

%% Load test set data PESHMELBA 
% the first simulation under each rain is considered the true yobs
cd(data_folder)
lhs_x_test = readNPY('LHS25_100_R500_dim6_sample_final.npy');
lhs_y_test = readNPY('Ymoisture_profile_LHS25_100_R500_dim6_YzeronS06_LogN.npy');

if (size(lhs_x_test,1) == size(lhs_y_test,1))
    n_lhs_test = size(lhs_x_test,1);
end

% take only the parameters of interest
expdes_x_test  = lhs_x_test(1:n_lhs_test, vector_input_indices);
expdes_x_train = lhs_x_train(1:n_lhs_train, vector_input_indices);

%% Calculate the cost function
depths = [0.5 1 2 3 4 5 6 10 15 20 25 30 35 40 45 50 55 65 75 100 150 200 250 300 400];
weighted = 0;
weights = ones(1,length(depths)); % depths

%% Reading input parameter marginal distributions, same for all parcels
cd(data_folder)
parcel_names = [545,530,480,503,526,527,481,524,525,529,544,539,528,523];
plot_idx = 4;
load('input_descriptor145.mat')

jj =1;
for ii = vector_input_indices
    InputOpts1.Marginals(jj).Type = deblank(Marginals.Type(ii,:));
    InputOpts1.Marginals(jj).Parameters = Marginals.Parameters{:,ii};
    jj = jj + 1;
end
myInput1 = uq_createInput(InputOpts1);

N_tot_train = 50; % the size of the experimental design with each rain in the train expdes

%% Train all the PCE metamodels
% if using R198, it is 54, if using test set R500, it is 1.
true_rain_idx = 1; % the observation comes from the "good rain", no forcing input error present, index starts from 1
true_test_idx = 1 + (true_rain_idx - 1) * N_tot_test; % index of true simulation, as present in test set
%%
for N_train = 50 %45:5:50 % number of simulation for PCE fitting
    for r = 1:500 %1:198 % number of cost functions to minimize

        % add Jb on the cost function !!!!!!!!!!
        expdes_J_train = ((lhs_y_train(1:n_lhs_train,:) - ones(n_lhs_train,1)*lhs_y_test(true_test_idx,:)).^2)*weights' + ...
            0.1 * 1/Marginals.Parameters{:,68}(2) *  (lhs_x_train(1:n_lhs_train,68)- Marginals.Parameters{:,68}(1)).^2+ ...
            0.1 * 1/Marginals.Parameters{:,71}(2) *  (lhs_x_train(1:n_lhs_train,71)- Marginals.Parameters{:,71}(1)).^2;

            name_file = strcat('pesh_profmoist_Jpce_errorLogN_Jb_Ridx', int2str(r));
            cd(code_folder)
            calculate_LARS_PCE_parcel(parcel_names(plot_idx),name_file, ...
                expdes_x_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train),:), ...
                expdes_J_train(((r-1)*N_tot_train+1):((r-1)*N_tot_train+N_train)),14,1, ...
                output_folder);
    end
end
% fix under which rain is the true yobs observed
true_test_theta = lhs_x_test(true_test_idx,vector_input_indices);
true_test_theta % check if its the true param. that was chosen % 0.3121 0.1836 0.3416 0.2827 0.0796 -3.6500

%% Get R2 on metamodels (validate on the test set available)
all_Q2_rains = [];
figure
for N_train = 50 % 45:5:50
    coeff_Q2_rains = [];
    R_temp = 200; % delete when all simulations available
    R = R_temp;
    for rain_idx = 1:500  
        % adapt cost function with Jb
        expdes_J_test = ((lhs_y_test(1:n_lhs_test,:) - ones(n_lhs_test,1)*lhs_y_test(true_test_idx,:)).^2)*weights'+ ...
            0.1 * 1/Marginals.Parameters{:,68}(2) *  (lhs_x_test(1:n_lhs_test,68)- Marginals.Parameters{:,68}(1)).^2+ ...
            0.1 * 1/Marginals.Parameters{:,71}(2) *  (lhs_x_test(1:n_lhs_test,71)- Marginals.Parameters{:,71}(1)).^2;
        
        cd(output_folder)
        load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN_Jb_Ridx', int2str(rain_idx),'_503'));
        Xval = expdes_x_test((2+(rain_idx-1)*N_tot_test):(rain_idx*N_tot_test),:); 
        Jval = expdes_J_test((2+(rain_idx-1)*N_tot_test):(rain_idx*N_tot_test));
        JLARS = uq_evalModel(myPCE_LARS,Xval);
        %figure
        hold on
        plot(Jval,JLARS,'bo')
        SSE = sum((Jval - JLARS).^2); 
        SST = sum((Jval-mean(Jval)).^2);
        coeff_R2 = 1 - SSE/SST;
        coeff_Q2_rains = [coeff_Q2_rains,coeff_R2];
    end
    all_Q2_rains = [all_Q2_rains, coeff_Q2_rains'];
end
plot(xlim,ylim,'-k', 'LineWidth',1)
%% Plot the Q2
%for i = 1:198
    hist(all_Q2_rains);
    find(all_Q2_rains < 0.90 ) % same bad metamodels as in 1st case, probably some bad simulations.
%end
%xlabel('N_{train}'); ylabel('R2'); title('R2 of each metamodel conditioned on false rain')
cd(output_folder)
writematrix(all_Q2_rains,'Q2_LARS_pesh_profmoist_Jpce_errorLogN_Jb_Ntrain_rep5_50_dim6_moistureprofile_YzeronS06_hourly.csv')

title('Q2 of the PCE metamodels on Ntrain = 50')
