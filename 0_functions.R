########################################################################
####       FUNCTIONS NOT DEPENDING ON ANY GLOBAL VARIABLES !!!!

#########################    PCA    ###############################
## transform GMM simulations to PCE coeffs
pca.recon <- function(pca, PC_scores, n_pc){
  mu <- matrix(rep(pca$center, nrow(PC_scores)), nrow = nrow(PC_scores), byrow = T)
  recon <- PC_scores[,1:n_pc] %*% t(pca$rotation[,1:n_pc]) + mu
  return(recon)
}

###########################     PCE EVALUATION     ###########################
isoprobTransfGaussian <- function(x,mean,stddev){
  return((x-mean)/stddev)
}

isoprobTransfGaussian_z_to_x <- function(z,mean,stddev){
  return(z*stddev+mean)
}

isoprobTransfUnif <- function(x,a,b){
  # transform U(a,b) ti U(-1,1)
  return(2*(x-a)/(b-a) - 1)
}

my.legendre <- function(x,n){
  # evaluates legendre polynomials of first n degrees in point x
  if (n==0){y <- 1}
  else if  (n==1){y <- x}
  else if  (n==2){y <- 1/2*(3*x^2 - 1)}
  else if  (n==3){y <- 1/2*(5*x^3 - 3*x)}
  else if  (n==4){y <- 1/8*(35*x^4 - 30*x^2 + 3)}
  else {stop("Error : Degree not supported")}
  return(y)}

## Fits the UQLab uqeval function !!
# eval_surrogate_OLS_in_point_x_pstcd(x, df_prior_Jb, A_multidx, coeff_PCE_fixed_rain) 
eval_surrogate_OLS_in_point_x_pstcd <- function(x, df_prior_Jb, A_multidx, coeff_PCE_fixed_rain) {
  # Isoprobabilistic transform of the point x
  unif_vec_idxs <- c(1,2,3,6,7); gauss_vec_idxs <-c(4,5,8) # dirty fix distrib. types
  # Lognormal variables to gaussians 4,8 : 
  x[4] <- log(x[4]); x[8] <- log(x[8])  
  Zgauss <- mapply(isoprobTransfGaussian, x[gauss_vec_idxs], df_prior_Jb$param1[gauss_vec_idxs], df_prior_Jb$param2[gauss_vec_idxs])
  Zunif <- mapply(isoprobTransfUnif, x[unif_vec_idxs], df_prior_Jb$param1[unif_vec_idxs], df_prior_Jb$param2[unif_vec_idxs])
  
  Psi_alpha = c()
  for (idx_alpha in seq(1,dim(A_multidx)[1])){
    alpha <- A_multidx[idx_alpha, ]
    alpha_gauss <- alpha[gauss_vec_idxs]
    alpha_unif <- alpha[unif_vec_idxs]
    
    # Gaussians ; hermite
    psi_alpha_Zgauss <- mapply(hermite, Zgauss, alpha_gauss, prob = TRUE)
    psi_alpha_Zgauss[alpha_gauss == 2] <- psi_alpha_Zgauss[alpha_gauss == 2]/sqrt(2) # normalize dirty
    psi_alpha_Zgauss[alpha_gauss == 3] <- psi_alpha_Zgauss[alpha_gauss == 3]/sqrt(2*3) # normalize dirty
    psi_alpha_Zgauss[alpha_gauss == 4] <- psi_alpha_Zgauss[alpha_gauss == 4]/sqrt(2*3*4) # normalize dirty
    
    # Uniforms; legendre
    psi_alpha_Zunif <- mapply(my.legendre, Zunif, alpha_unif)
    psi_alpha_Zunif[alpha_unif == 1] <- psi_alpha_Zunif[alpha_unif == 1]*sqrt(2*1+1) # normalize dirty
    psi_alpha_Zunif[alpha_unif == 2] <- psi_alpha_Zunif[alpha_unif == 2]*sqrt(2*2+1) # normalize dirty
    psi_alpha_Zunif[alpha_unif == 3] <- psi_alpha_Zunif[alpha_unif == 3]*sqrt(2*3+1) # normalize dirty
    psi_alpha_Zunif[alpha_unif == 4] <- psi_alpha_Zunif[alpha_unif == 4]*sqrt(2*4+1) # normalize dirty
    
    # Combine all in one vector psi_alpha_Z
    psi_alpha_Z <- c()
    psi_alpha_Z[unif_vec_idxs] <- psi_alpha_Zunif; 
    psi_alpha_Z[gauss_vec_idxs] <- psi_alpha_Zgauss;
    psi_alpha_Z
    
    Psi_alpha_Z <- prod(psi_alpha_Z)
    Psi_alpha <- cbind(Psi_alpha, Psi_alpha_Z)
    Psi_alpha
  }    
  colnames(Psi_alpha) <- paste(A_multidx[seq(1,dim(A_multidx)[1]),1], 
                               A_multidx[seq(1,dim(A_multidx)[1]),2],
                               A_multidx[seq(1,dim(A_multidx)[1]),3], 
                               A_multidx[seq(1,dim(A_multidx)[1]),4],
                               A_multidx[seq(1,dim(A_multidx)[1]),5],
                               A_multidx[seq(1,dim(A_multidx)[1]),6],
                               A_multidx[seq(1,dim(A_multidx)[1]),7],
                               A_multidx[seq(1,dim(A_multidx)[1]),8])
  # Psi_alpha is a row vector, coeff_PCE_fixed_rain is a column vector
  J_surr_OLS <- Psi_alpha%*%coeff_PCE_fixed_rain
  return(J_surr_OLS)
}

eval_surrogate_OLS_in_point_x_variable_rain_pstcd <- function(x, df_prior, A_multidx, ML_coeff_multidx_r)
{
  values <- c()
  for (rain_idx in seq(1,dim(ML_coeff_multidx_r)[2]))
  {
    #print(paste('rain = ', rain_idx))
    coeff_PCE_fixed_rain <- ML_coeff_multidx_r[,rain_idx]
    values <- c(values,  eval_surrogate_OLS_in_point_x_pstcd(x, df_prior, A_multidx, coeff_PCE_fixed_rain))
  }
  return(values)
}

eval_surrogate_OLS_in_point_x <- function(x, df_prior, A_multidx, coeff_PCE_fixed_rain) {
  # evaluates the PCE in a fixed 6-dimensional point x
  Z <- mapply(isoprobTransfGaussian, x, df_prior$mean, df_prior$sd)
  #print('evaluate OLS normalized hermite')
  Psi_alpha = c()
  for (idx_alpha in seq(1,dim(A_multidx)[1])){
    alpha <- A_multidx[idx_alpha, ]
    psi_alpha_Z <- mapply(hermite, Z, alpha, prob = TRUE)
    psi_alpha_Z[alpha == 2] <- psi_alpha_Z[alpha == 2]/sqrt(2) # normalize dirty
    psi_alpha_Z[alpha == 3] <- psi_alpha_Z[alpha == 3]/sqrt(2*3) # normalize dirty
    psi_alpha_Z[alpha == 4] <- psi_alpha_Z[alpha == 4]/sqrt(2*3*4) # normalize dirty
    
    Psi_alpha_Z <- prod(psi_alpha_Z)
    Psi_alpha <- cbind(Psi_alpha, Psi_alpha_Z)
  }    
  # exp des x card of basis
  colnames(Psi_alpha) <- paste(A_multidx[seq(1,dim(A_multidx)[1]),1], 
                               A_multidx[seq(1,dim(A_multidx)[1]),2],
                               A_multidx[seq(1,dim(A_multidx)[1]),3], 
                               A_multidx[seq(1,dim(A_multidx)[1]),4],
                               A_multidx[seq(1,dim(A_multidx)[1]),5],
                               A_multidx[seq(1,dim(A_multidx)[1]),6])
  
  J_surr_OLS <- Psi_alpha%*%coeff_PCE_fixed_rain
  return(J_surr_OLS)
}

eval_surrogate_OLS_in_point_x_variable_rain <- function(x, df_prior, A_multidx, ML_coeff_multidx_r)
{
  values <- c()
  for (rain_idx in seq(1,dim(ML_coeff_multidx_r)[2]))
  {
    #print(paste('rain = ', rain_idx))
    coeff_PCE_fixed_rain <- ML_coeff_multidx_r[,rain_idx]
    values <- c(values,  eval_surrogate_OLS_in_point_x(x, df_prior, A_multidx, coeff_PCE_fixed_rain))
  }
  return(values)
}

###########################   READING FILES      ###########################
## Read .npy peshmelba simulations to dataframe of cost function, unweighted
read_npy_peshmelba_to_dataframe_J_squared_all_factors <- function (input_sample_final_train,Ymoisture_profile_train,
                                                                   input_sample_final_test,Ymoisture_profile_test,
                                                                   true_test_idx, folder)
{
  setwd(peshmelba_simulations_folder)
  Y_test  <- npyLoad(Ymoisture_profile_test)
  Y_true <- Y_test[true_test_idx,]
  setwd(folder)
  Y_train <- npyLoad(Ymoisture_profile_train)
  x_train <- npyLoad(input_sample_final_train)
  if (length(Y_train[1,]) != length(Y_true))
  {stop("Error : dimensions not matching")}
  delta_Y <- Y_train - t(matrix(Y_true, length(Y_true),length(Y_train[,1])))
  delta_Y2 <- delta_Y^2
  weights <- rep(1,length(Y_train[1,]))
  J_train <- delta_Y2%*%weights
  x_J_dataframe <- data.frame(x_train, J = J_train)
  return(x_J_dataframe)
}

## Read .npy peshmelba simulations to dataframe of cost function, unweighted
read_npy_peshmelba_to_dataframe_J_squared_6factors_folder <- function (input_sample_final_train,Ymoisture_profile_train,
                                                                       input_sample_final_test,Ymoisture_profile_test,
                                                                       true_test_idx, folder)
{
  setwd(folder)
  Y_test  <- npyLoad(Ymoisture_profile_test)
  Y_true <- Y_test[true_test_idx,]
  Y_train <- npyLoad(Ymoisture_profile_train)
  x_train <- npyLoad(input_sample_final_train)
  if (length(Y_train[1,]) != length(Y_true))
  {stop("Error : dimensions not matching")}
  delta_Y <- Y_train - t(matrix(Y_true, length(Y_true),length(Y_train[,1])))
  delta_Y2 <- delta_Y^2
  #  weights = c(0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300,400)
  weights <- rep(1,length(Y_train[1,]))
  J_train <- delta_Y2%*%weights
  x_J_dataframe <- data.frame(th9   = x_train[,63], 
                              mn10  = x_train[,68], 
                              th10  = x_train[,72], 
                              th13  = x_train[,99], 
                              thr10 = x_train[,71], 
                              hg10  = x_train[,65], 
                              J = J_train) # ML order
  return(x_J_dataframe)
}

## Read .npy peshmelba simulations to dataframe of cost function, WEIGHTED
read_npy_peshmelba_to_dataframe_J_squared_6factors_folder_weighted <- function (input_sample_final_train,Ymoisture_profile_train,
                                                                                input_sample_final_test,Ymoisture_profile_test,
                                                                                true_test_idx, folder)
{
  setwd(folder)
  Y_test  <- npyLoad(Ymoisture_profile_test)
  Y_true <- Y_test[true_test_idx,]
  Y_train <- npyLoad(Ymoisture_profile_train)
  x_train <- npyLoad(input_sample_final_train)
  if (length(Y_train[1,]) != length(Y_true))
  {stop("Error : dimensions not matching")}
  delta_Y <- Y_train - t(matrix(Y_true, length(Y_true),length(Y_train[,1])))
  delta_Y2 <- delta_Y^2
  # depths of the upper and lower bound of each cell
  depth_upper = c(0,0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300)
  depth_lower = c(0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300,400)
  width = (depth_lower - depth_upper) # width of each cell
  weights <- width
  J_train <- delta_Y2%*%weights
  x_J_dataframe <- data.frame(th9   = x_train[,63], 
                              mn10  = x_train[,68], 
                              th10  = x_train[,72], 
                              th13  = x_train[,99], 
                              thr10 = x_train[,71], 
                              hg10  = x_train[,65], 
                              J = J_train) # ML order
  return(x_J_dataframe)
}

###########################   PLOTTING      ###########################
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, breaks = 10, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}

###########################   GGPLOT TO TIKZ

ggplot_tikzDevice <- function(ggplot, directory, name, width, height)
{
  library(ggplot2)
  library(grid)
  require(tikzDevice)
  
  setwd(directory)
  tikz(paste(name,'.tex', sep=""), standAlone = TRUE, width=width, height=height)
  #gg_arranged_traj_sobol_rain
  ggplot
  #Close the device
  dev.off()
  
  # Compile the tex file
  tools::texi2dvi(paste(name,'.tex', sep=""),pdf=T)
  
  # optionally view it:
  system(paste(getOption('pdfviewer'),paste(name,'.tex', sep="")))
}




