%% Get the basis union A
clearvars
rng(100,'twister')
cd '/home/katarina.radisic/UQLab_Rel2.0.0/core/'
uqlab
output_folder='/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels';
code_folder = '/home/katarina.radisic/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/2_Robustness_scripts_Abasis';
data_folder = '/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw';

%% Set data
R_train=200; % nuber of rains used to obtain the training capha

%% Get basis Ar
beta = 0.00;
str_beta = '00';

a_old = []; c_old = [];
N_train = 50; n_inputs = 6;
% try Ntrain = 40, 45, 50
for rain_idx = 1:R_train
cd(strcat(output_folder))
    %load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN_Jb_Ridx', int2str(rain_idx),'_503'))
    load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',int2str(beta),'_Ridx', int2str(rain_idx),'_503'))
    cd(code_folder)
    [a_new, c_new] = get_multiindex(myPCE_LARS,n_inputs);
    [a_old, c_old] = union_multiindices(a_old,c_old,a_new,c_new);
end

%% Reduce the basis
cd(code_folder)
[a_reduced, c_reduced] = discard_small_regressors(a_old,c_old,10^(-10),n_inputs);

%a_reduced = a_old;

%% save
cd(output_folder)
vector_input_indices = [63,68,72,99,71,65];

save(strcat('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.mat'),"a_reduced")

T = array2table(a_reduced);
T.Properties.VariableNames(1:6) = split(int2str(vector_input_indices))
writetable(T,strcat('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv'))
