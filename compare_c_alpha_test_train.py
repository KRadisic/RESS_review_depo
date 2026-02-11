import pandas as pd
import matplotlib.pyplot as plt


beta = 0.005
str_beta ='005'

ml_coeff_multidx_r_TRAIN = pd.read_csv(f"ml_multicoeff_c_alpha_OLS_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb{str_beta}.csv", header = None)
ml_coeff_multidx_r_TEST = pd.read_csv(f"ml_multicoeff_c_alpha_OLS_TEST_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb{str_beta}.csv", header = None)

#bins = numpy.linspace(-10, 10, 100)

for ii in range(17) :
    plt.hist(ml_coeff_multidx_r_TRAIN.loc[ii], alpha=0.5, label='x', density=True)
    plt.hist(ml_coeff_multidx_r_TEST.loc[ii], alpha=0.5, label='y', density=True)
    plt.legend(loc='upper right')
    plt.show()


 # 200 realizations of c_alpha
# 500 realizations of c_alpha
# 10 14 23 43 73 80 82 86 87 123 130 131 151 155 164 173 175 178 188

