%% Export PCE coefficients in csv format
clearvars
rng(100,'twister')
cd '/home/kradisic/Documents/my_UQLab_folder/UQLab_Rel2.2.0/core' 

uqlab
output_folder='/home/kradisic/Documents/RESS_review_depo/metamodels';
data_folder = '/home/kradisic/Documents/RESS_review_depo/';

%% Set data
R_train=500;
beta=0.005;
str_beta = '005';
N_train = 50;

%% Read TRAIN metamodels OLS
c_alpha = [];
for rain_idx = 1:200
    cd(output_folder)

    load(strcat('myPCE_OLS_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'_Ridx',int2str(rain_idx), 'N_train', int2str(N_train), '_503.mat'),'myPCE_OLS');
    c_alpha = [c_alpha, myPCE_OLS.PCE.Coefficients];
end

%% pair plot 
%figure
%plotmatrix((c_alpha(:,:))')

%% Save the c_alpha
cd(output_folder)
writematrix(c_alpha,strcat('ml_multicoeff_c_alpha_OLS_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv')) 

%% Read VALIDATION metamodels OLS
c_alpha_validation = [];
for rain_idx = 1:499
    cd(output_folder)
    load(strcat('myPCE_OLS_TEST_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'_Ridx',int2str(rain_idx),'_503.mat'),'myPCE_OLS');
    c_alpha_validation = [c_alpha_validation, myPCE_OLS.PCE.Coefficients];
end

%% pair plot
%figure
%plotmatrix((c_alpha_validation(:,:))')

%% Save the c_alpha
cd(output_folder)
writematrix(c_alpha_validation,strcat('ml_multicoeff_c_alpha_OLS_TEST_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv')) 