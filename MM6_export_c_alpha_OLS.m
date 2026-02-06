%% 
clearvars
rng(100,'twister')
cd '/home/katarina.radisic/UQLab_Rel2.0.0/core/'
uqlab
output_folder='/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels';
code_folder = '/home/katarina.radisic/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/2_Robustness_scripts_Abasis';
data_folder = '/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw';

%% Set data
R_train=500;
beta=0.00;
str_beta = '00';

%% Read TRAIN metamodels OLS
c_alpha = [];
for rain_idx = 1:200
    cd(output_folder)
    %load(strcat('myPCE_OLS_pesh_profmoist_Jpce_errorLogN_Jb_Ridx',int2str(rain_idx),...
    %        '_503.mat'),'myPCE_OLS');
    load(strcat('myPCE_OLS_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'_Ridx',int2str(rain_idx),'_503.mat'),'myPCE_OLS');
    c_alpha = [c_alpha, myPCE_OLS.PCE.Coefficients];
end

%% pair plot
figure
plotmatrix((c_alpha(:,:))')

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
figure
plotmatrix((c_alpha_validation(:,:))')

%% Save the c_alpha
cd(output_folder)
writematrix(c_alpha_validation,strcat('ml_multicoeff_c_alpha_OLS_TEST_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv')) 