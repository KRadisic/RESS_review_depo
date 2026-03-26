####################
## MATLAB scripts (using UQLab)
####################

MM1_figure_quality_PCE_Q2.m
MM1_updated_Robustness_script_R2_external_true.m
Train Least Angle Regression (LARS) Polynomial Chaos Expansion (PCE) metamodels and obtain their Q2 on the test set. Compare for increasing number of training points. Save the Q2 data for pretty plotting in R the Figure 6

MM2_Get_union_A_basis.m
Runs through all bases obtained with LARS, obtains their union, removes PCE components with small coefficients (controlled by a chosen threshold).
 
MM3_Get_OLS_train_validation.m
Fit Ordinary Least Squares (LARS) Polynomial Chaos Expansion (PCE) metamodels on the basis fixed in the previous step, for both the train experimental design and trajectories, and the test experimental design and trajectories. Save all PCE metamodels in UQLab format (the c_alpha coefficients for the train and test set are exported in MM7)

MM4_Get_LARS_minima.m
Minimization with BFGS algirithm of individual train PCE trajectories.

MM5_export_c_alpha_OLS.m
Opens the OLS PCE metamodels obtained in MM3, exports and saves their coefficients (c_alpha train and test) in .csv format


####################
## R scripts 
####################

7_get_PCA_from_train_set_sigma02.R
Read the PCE cefficients c_alpha exported from matlab, perform principal component analysis (PCA) on the c_alpha from the train set, infer the distribution of the projected coefficients with kernel density estimates, generate nes samples c_alpha new, read the c_alpha from test set exported from matlab. sript used for Figure 7 and Figure 8. 

8_plot_eps_marg_valid_beta.R
Get the eps_marg obtained through bootstrapping of the validation set, then used for the gray area on convergence figures for eps_marg

9_plot_eps_marg_Ntrain_boot.R
Script used to evaluate the eps_marg error of the stochastic emulator for increasing number of Ntrain points, script used for Figure 15.

9_plot_eps_marg_Rtrain_boot.R
Script used to evaluate the eps_marg error of the stochastic emulator for increasing number of Rtrain trajectories, script used for Figure 9.

10_plot_new_cost_function_and_profiles.R
Read the new moisture profiles obtained with PESHMELBA simulations with the robust estimators, compare with classic calibration through the cost function and visually. Script used for Figure 14.

11_plot_convergence_robust_estim_nLHS.R
Implements Algorithm 1, which evaluates the stochastic emulator in an increasing number of candidate points. This is done for an increasing number of simulations of the stochastic emulator. (a parallelized version of the script for use on cluster is 11_plot_convergence_robust_estim_nLHS_parallelized.R)

11_plot_convergence_robust_estim_nLHS_thresholds_posttraitement.R
The robust estimators are obtained starting from the evaluations of the stochastic emulator obtained in 11_plot_convergence_robust_estim_nLHS.R, these  are used to obtain robust estimators for an increasing number of candidate points and for different thresholds. Script used for Figure 12 illustrating the probability of exceeding the fixed threshold with respect to the chosen threshold.

12_plot_robust_estim_histograms_make_PESHMELBA_LHS.R
Arrange results from robust calibration in the input format needed for PESHMELBA simulations. Plot the results of classical calibration and robust calibration, used for Figure 11.

13_plot_MM_stoch_realizations.R
Evaluate the stochastic metamodel in the neighbourhood of the robust estimator, script used to obtain Figure 13. Arrange and save results from robust calibration in the input format needed for PESHMELBA simulations.

13_temp_plot_MM_stoch_vs_PESHrealizations.R
Compare visually trajectories obtained from PESHMELBA simulations to realizations of the stochastic emulator in the neighbourhood of the robust estimator.

illustration_KDE.R
Script for annex example for Kernel Density Estimate

visualize_calpha_boot.R
Script used only for Figure A4 in the letter to the reviewer, 

visualize_raw_input_data_RESS.R
Visualizes the scatterplot of raw input datasets of the train and test set. Script used only for Figure A1 in the letter to the reviewer.


####################
## MATLAB functions
####################

calculate_LARS_PCE_parcel.m
Evaluates and saves a PCE metamodel with Least Angle Regression (calls uq_createModel function from UQLab)

calculate_OLS_PCE_parcel.m
Evaluates and saves a PCE metamodel on a fixed basis A with Ordinary Least Squares (calls uq_createModel function from UQLab)

discard_small_regressors.m
Remove from the basis the PCE components whose corresponding coefficients are smaller than a given threshold eps.

get_multiindex.m
Get the indices from a PCE metamodel object

readNPY.m
readNPYheader.m
Both needed for reading .npy Python data in matlab

union_multiindices.m
Get the union of multiple PCE bases.


#########################
###   DATA 
#########################

input_descriptor145.mat
Distributions of the 145 PESHMELBA parameters in the format used by UQLab

entry_factor_names.csv
Names of the 145 PESHMELBA parameters. Only input_indices  63, 65, 68, 71, 72, 99 are variable in this work.

prior_params_wo_bds145.csv
(Hyper)parameters of the 145 PESHMELBA parameter distributions in csv format. Only input_indices  63, 65, 68, 71, 72, 99 are variable in this work.

union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb005.mat 
union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb005.csv
The fixed PCE basis used for the Ordinary Least Squares PCE (17 components) both in Matlab format from UQlab and in csv.

union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb005_boot.csv
Union of all the PCE bases obtained with Least Angle Regression PCE (32 components), before removing the component with small regressors.

