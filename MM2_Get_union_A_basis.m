%% Get the basis union A
clearvars
rng(100,'twister')
cd /Users/icj_ecl/Documents/UQLab_Rel2.1.0/core/
uqlab
output_folder='/Users/icj_ecl/Documents/RESS_review_depo/metamodels/';
data_folder = '/Users/icj_ecl/Documents/RESS_review_depo/'; %/data/raw

%% Set data
R_train=200; % nuber of rains used to obtain the training capha

%% Get basis Ar
beta = 0.005;
str_beta = '005';

a_old = []; c_old = [];
N_train = 50; n_inputs = 6;
% try Ntrain = 40, 45, 50
for rain_idx = 1:R_train
cd(output_folder)
    %load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN_Jb_Ridx', int2str(rain_idx),'_503'))
    load(strcat('myPCE_LARS_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb',str_beta,'_Ridx', int2str(rain_idx),'_503'))
    cd(data_folder)
    [a_new, c_new] = get_multiindex(myPCE_LARS,n_inputs);
    [a_old, c_old] = union_multiindices(a_old,c_old,a_new,c_new);
end

%% Reduce the basis
cd(data_folder)
[a_reduced, c_reduced] = discard_small_regressors(a_old,c_old,10^(-10),n_inputs);

a_reduced = a_old;

%% save
cd(output_folder)
vector_input_indices = [63,68,72,99,71,65];

save(strcat('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.mat'),"a_reduced")

T = array2table(a_reduced);
T.Properties.VariableNames(1:6) = split(int2str(vector_input_indices))
writetable(T,strcat('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv'))
