%% METAMODEL AND MINIMIZE THE MEAN OF THE COST FUNCTIONS WITH PCE
clearvars
rng(100,'twister')
cd '/home/katarina.radisic/UQLab_Rel2.0.0/core/'
uqlab
output_folder='/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels';
code_folder = '/home/katarina.radisic/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/2_Robustness_scripts_Abasis';
data_folder = '/home/katarina.radisic/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw';

%% Set data
vector_input_indices = [63,68,72,99,71,65];
N_tot_test = 100; % number of simulations after which the rain changes 
N_tot_train = 50; % the size of the experimental design with each rain in the train expdes


%% Minimisation of the mean with the PCE LARS metamodel, do we find theta star ?

cd(output_folder)
load('myPCE_LARS_pesh_profmoist_Jpce_errorLogN_Jb_MEAN_Ridx200_503');
evaluate_pce_handle = @(x) uq_evalModel(myPCE_LARS,x);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.229577, 0.1237539, 0.2183556, 0.2332047, 0.01391467, -3.6]; % fix quantiles
ub = [0.434923, 0.2344461, 0.4136644, 0.4417953, 0.10844333, -2]; % fix quantiles
x0 = [0.31,0.20,0.34,0.28,0.07,-3.6];
[xmin,fval,exitflag,output] = fmincon(evaluate_pce_handle,x0,A,b,Aeq,beq,lb,ub);  
X_min_LARS = xmin
f_min_LARS = fval

%% save the minimum of the mean
cd(output_folder)
vector_input_names = ["thetas9","mn10","thetas10","thetas13","thetar10","hg10"];
X_min_mean = [[vector_input_names,"fcondmin"]; [xmin,fval]];

writematrix(X_min_mean,strcat('mean_minimizers_LARS_errorLogN_Jb_',name_train,'_Ntrain_hourly_BDDquantiles',...
    int2str(N_train),'.csv'));
