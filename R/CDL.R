#' @title Compressive Dictionary Learning
#' @description Implementation of the Sketched Stochastic Dictionary Learning (SSDL) method, which 
#' learns a dictionary from a large-scale matrix-like dataset by operating on a compressed version of data (a.k.a. data sketch). 
#' @details 
#' \code{CDL} builds a dictionary by alternating two steps: 
#' calculating the code matrix that contains the weights of the dictionary elements, 
#' and updating the dictionary. For computational efficiency, the code matrix is 
#' computed only for a randomly selected subset of data vectors \eqn{x_1, \dots, x_n} 
#' (a.k.a. batch). The code matrix is obtained as a solution of the following optimization problem:
#' \eqn{\min\limits_{A \in R^{+}_{K\times n}} \sum_{i=1}^{n} \|x_i - D\cdot \alpha_i\|^2 + 
#' \lambda \cdot\|\alpha_i\|_1}, where \eqn{n} denotes a batch size,
#' \eqn{A = \{\alpha_1,\dots, \alpha_{n}\}} is a code matrix and \eqn{\lambda} is a regularization parameter which
#' defines the data sparsity level.
#' 
#' The dictionary is updated by taking one step along the gradient of the objective function
#' \eqn{F(D, A) = \|SK(Data) - SK(A\cdot D)\|^2}. 
#' Two gradient descent update rules are available: Nesterov accelerated and momentum. 
#' 
#' \eqn{SK(\cdot)} is a sketch operator, which compresses a matrix 
#' into a fixed size complex vector referred to as a data sketch.  It has been introduced in \insertRef{DBLP:journals/corr/KerivenBGP16}{SSDL} and 
#' it is defined as \eqn{SK(Data) = \frac{1}{N}\sum_{i=1}^N \exp{\left(-1i \cdot W\cdot x_i\right)}},
#'  where \eqn{W} is a frequency matrix and \eqn{x_1, \dots, x_N} are data vectors.   
#' The data compression is performed using routines from [chickn](https://CRAN.R-project.org/package=chickn) package. 
#' 
#' \code{CDL} allows also to use the decaying learning rate, \emph{i.e.} 
#' \eqn{\code{learn\_rate}^t = \frac{\code{learn\_rate}}{1+ (t-1)\cdot \code{gamma}}}, where \eqn{t} is the iteration number. 
#' 
#' @param Data is a Filebacked Big Matrix \eqn{s \times N} with data vectors stored in the matrix columns.
#' @param K is a dictionary size. 
#' @param D is an initial dictionary. If it is NULL, the dictionary is initialized by
#' random selection of \code{K} signals from \code{Data} and it is saved in the \code{DIR_tmp} directory.
#' @param SK_Data is a data sketch. It is a \eqn{2m}-dimensional complex vector. The first \eqn{m} coordinates 
#' correspond to the real parts and the last \eqn{m} coordinates to the imaginary parts.  
#' If it is NULL, the sketch is computed using \code{Sketch} function of 
#' [chickn](https://CRAN.R-project.org/package=chickn) package.
#' @param Frequencies is a frequency matrix \eqn{m \times s} with frequency vectors in the matrix rows. 
#' If NULL, the frequencies are generated using \code{GenerateFrequencies} 
#' function of [chickn](https://CRAN.R-project.org/package=chickn) package. 
#' @param m is a number of the frequency vectors.
#' @param pos.dic indicates whether the dictionary is positive (default) or not. 
#' @param maxEpoch is a number of epochs. 
#' @param batch_size is a batch size.
#' @param lambda is a regularization parameter.
#' @param typeOptim is a type of the optimization scheme used in the dictionary update. 
#' Possible values are c('Nesterov', 'Momentum'). It is suggested to use 'Nesterov' scheme.
#' @param learn_rate is a learning rate value. The default value is 0.1.
#' @param gamma is a decay parameter. The default value is 0, which corresponds to the constant learning rate.
#' @param alpha is a momentum weight.
#' @param ncores is a number of cores. The default value is 1. 
#' @param DIR_tmp is a directory to save the initial dictionary and intermediate results.
#' @param grad_t_1 is an initial momentum matrix. By default it is NULL, and it is initialized as a zero matrix. 
#' @param verbose controls how much output is shown and saved during the optimization process. Possible values:
#' \itemize{
#'     \item 0 no output (default value)
#'     \item 1 show iteration number and value of objective function
#'     \item 2 1 + save a dictionary and a momentum matrix at the end of each epoch.
#' }
#' @param ... are additional parameters passed to \link[chickn]{GenerateFrequencies} function.
#' @return a list
#'     \itemize{
#'        \item \code{D} is the obtained dictionary,
#'        \item  \code{objFunProcess} is objective function values computed at the end of each iteration,
#'        \item \code{learning_rate} is learning rate values.
#'     }
#' @examples 
#' \donttest{
#' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
#' X_fbm = bigstatsr::as_FBM(X)$save()
#' W = chickn::GenerateFrequencies(Data = X_fbm, m = 64, N0 = ncol(X_fbm),
#'                                 ncores = 1, niter= 3, nblocks = 2, sigma_start = 0.001)$W
#' SK= chickn::Sketch(X_fbm, W)
#' D = CDL(Data = X_fbm, K = 10, SK_Data = SK, Frequencies = W,
#'         D = NULL, pos.dic = TRUE, maxEpoch = 3, batch_size = 50, 
#'         lambda = 0, learn_rate = 0.1, alpha = 0.9, 
#'         gamma = 0, ncores = 2, DIR_tmp = tempdir(), 
#'         verbose=0, typeOptim = "Nesterov")$D
#'}                   
#' @importFrom chickn Sketch GenerateFrequencies
#' @importFrom parallel mclapply makeCluster stopCluster 
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#' @references 
#' \itemize{
#'       \item \insertRef{permiakova2021SSDL}{SSDL}
#'       \item \insertRef{permiakova2021chickn}{SSDL}}
#' @seealso \code{\link{Gradient_D_cpp_parallel}}, [chickn](https://CRAN.R-project.org/package=chickn), 
#' \code{chickn::Sketch}, \code{chickn::GenerateFrequencies}
#' 
#' @export
CDL <-function(Data, K,
               SK_Data = NULL, 
               Frequencies = NULL, 
               D = NULL, pos.dic = TRUE, 
               learn_rate = 0.1, alpha = 0.9, gamma = 0, 
               maxEpoch = 5, batch_size, lambda = 0,
               ncores = 1,
               typeOptim = 'Nesterov', 
               DIR_tmp = tempdir(), grad_t_1 = NULL,
               verbose = 0, m = nrow(Frequencies), ...){
  
  # ========= Initialization ================================
  RcppParallel::setThreadOptions(numThreads = ncores)
  d = nrow(Data)
  N = ncol(Data)
  nIter = maxEpoch*ceiling(N/batch_size)
  if(length(Frequencies)==0){
    Frequencies = chickn::GenerateFrequencies(Data = Data, m = m, ...)$W
  }
  if(length(SK_Data) ==0){
    SK_Data = chickn::Sketch(Data = Data, W = Frequencies, ind.col = 1:ncol(Data), ncores = ncores, parallel = FALSE)    
  }
  #m = nrow(Frequencies)
  if(length(D)==0){
    D = matrix(NA, ncol = K, nrow = d)
    D <- Data[, sample(ncol(Data), K)]
    save(D, file = file.path(DIR_tmp, 'InitDict.RData'))
    t0=1
  }else{
    if(ncol(D)== K){
      t0=1
    }else{
      stop(paste0('Error: incorrect initial dictionary dimension. Require ', K, ' dictionary atoms'))
    }
  }
  if(verbose > 0){
    value_OF = rep(NA, nIter)
#    value_OF_init = rep(NA, nIter) 
  }
  grad_D = matrix(NA, ncol = K, nrow = d)
  if(length(grad_t_1)==0){
    grad_t_1= matrix(0, ncol = K, nrow = d)
  }
  
  A = matrix(NA, ncol = batch_size, nrow = K)#, sparse = TRUE)

  lr = learn_rate*(1/(1+(0:nIter)*gamma))
  alpha = rep(alpha, nIter)
  
  if(pos.dic){
    lb_dic = rep(0, nrow(D))
  }else{
    lb_dic = rep(-Inf, nrow(D))
  }
  #====================================== MAIN LOOP =============================================================    
  switch(typeOptim, 
         'Nesterov' = {
           for(epoch in 1:maxEpoch){
             ind_sample = matrix(sample(N, N, replace = FALSE), nrow = batch_size) 
             for(t in 1:ncol(ind_sample)){
               if(verbose == 1){cat("Iteration", t + t0-1, '\n')}
               if(Sys.info()['sysname'] == 'Linux' || Sys.info()['sysname'] == 'SunOS'){
                 A = simplify2array(parallel::mclapply(ind_sample[,t],
                                                       function(i){
                                                         as.vector(glmnet::glmnet(x = D, y = Data[,i],
                                                                                  lambda = lambda,
                                                                                  intercept = FALSE, lower.limits = 0,
                                                                                  standardize = T)$beta)
                                                       }, mc.cores = ncores))
               }else{
                if(Sys.info()['sysname'] == 'Windows' || Sys.info()['sysname'] == 'Darwin'){
                  core <- parallel::makeCluster(ncores, type = 'PSOCK')
                  doParallel::registerDoParallel(core)
                  doRNG::registerDoRNG()
                  on.exit(parallel::stopCluster(core))
                  i<-NULL
                  A <-foreach::foreach(i = ind_sample[,t], .combine = 'cbind') %dopar% {
                    as.vector(glmnet::glmnet(x = D, y = Data[,i],
                                             lambda = lambda,
                                             intercept = FALSE, lower.limits = 0,
                                             standardize = T)$beta)
                  }
                }else{
                  stop(paste("Package is not compatible with", Sys.info()['sysname']))
                } 
               }
               #Dictionary updata
               # if(OF_compute){
               #   output = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data, ComputeGrad = F)
               #   value_OF_init[t0+t-1] = output$ObjFun
               # }
               
               D = D - alpha[t+t0-1]*grad_t_1
               D = apply(D, 2, function(x, lb) pmax(lb, x), lb = lb_dic)
               
               output = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data)
               grad_D = output$grad
               
               D = D - lr[t0+t-1]*grad_D
               grad_t_1 = alpha[t0+t-1] * grad_t_1 + lr[t0+t-1]*grad_D # 
               
               #grad_t_1 = alpha[t+t0-1] * grad_t_1 - lr[t0+t-1]*grad_D # CORRECT FORMULAS
               #out2 = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data, ComputeGrad = FALSE)$ObjFun
               D = apply(D, 2, function(x, lb) pmax(lb, x), lb = lb_dic)
               D = D/matrix(pmax(rep(1, ncol(D)), sqrt(colSums(D^2))), ncol = ncol(D), nrow = nrow(D), byrow = T)
               if(verbose > 0){
                 value_OF[t0+t-1] = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data, ComputeGrad = FALSE)$ObjFun
                 cat("Objective Function: ", value_OF[t0+t-1], '\n')
               }
               
               # if((t0+t-1)%% saveEach == 0){
               #   save(D, file =  file.path(DIR_tmp, paste0('Dictionary_iter', t0+t-1, '.RData')))
               # }
             }
             t0=t0+t
             if(verbose > 1){
               save(D, file =  file.path(DIR_tmp, paste0('Dictionary_epoch', epoch, '.RData')))
               save(grad_t_1, file =  file.path(DIR_tmp, paste0('Momentum_matrix_epoch', epoch, '.RData')))
               }
           }
         },
         'Momentum' = {
           for(epoch in 1:maxEpoch){
             ind_sample = matrix(sample(N, N, replace = FALSE), nrow = batch_size) 
             for(t in 1:ncol(ind_sample)){
               if(verbose > 0){cat("Iteration", t + t0-1, '\n')}
               if(Sys.info()['sysname'] == 'Linux' || Sys.info()['sysname'] == 'SunOS'){
                 A = simplify2array(parallel::mclapply(ind_sample[,t],
                                                       function(i){
                                                         as.vector(glmnet::glmnet(x = D, y = Data[,i],
                                                                                  lambda = lambda,
                                                                                  intercept = FALSE, lower.limits = 0,
                                                                                  standardize = T)$beta)
                                                       }, mc.cores = ncores))
               }else{
                 if(Sys.info()['sysname'] == 'Windows' || Sys.info()['sysname'] == 'Darwin'){
                   core <- parallel::makeCluster(ncores, type = 'PSOCK')
                   doParallel::registerDoParallel(core)
                   doRNG::registerDoRNG()
                   on.exit(parallel::stopCluster(core))
                   i<-NULL
                   A <-foreach::foreach(i = ind_sample[,t], .combine = 'cbind') %dopar% {
                     as.vector(glmnet::glmnet(x = D, y = Data[,i],
                                              lambda = lambda,
                                              intercept = FALSE, lower.limits = 0,
                                              standardize = T)$beta)
                   }
                 }else{
                   stop(paste("Package is not compatible with", Sys.info()['sysname']))
                 } 
               }
               # A = simplify2array(parallel::mclapply(ind_sample[,t],
               #                                       function(i){
               #                                         if(sum(Data[, i]) >0){
               #                                           as.vector(glmnet::glmnet(x =D, 
               #                                                                    y = Data[,i],
               #                                                                    lambda = lambda,
               #                                                                    intercept = FALSE, lower.limits = 0,
               #                                                                    standardize = T)$beta)
               #                                         }else{
               #                                           rep(0, ncol(D))
               #                                         }
               #                                       }, mc.cores = ncores))
               output = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data, ComputeGrad = TRUE)
               grad_D = output$grad
               grad_t_1 = alpha[t0+t-1] * grad_t_1 + lr[t0+t-1]* grad_D
               D = D -grad_t_1
               D = apply(D, 2, function(x, lb) pmax(lb, x), lb = lb_dic)
               D = D/matrix(pmax(rep(1, ncol(D)), sqrt(colSums(D^2))), ncol = ncol(D), nrow = nrow(D), byrow = T)
               
               if(verbose >0){
                 value_OF[t0+t-1] = Gradient_D_cpp_parallel(D, A, Frequencies, SK_Data, ComputeGrad = FALSE)$ObjFun
                 cat("Objective Function: ", value_OF[t0+t-1], '\n')
                }
               # if((t0+t-1)%% saveEach == 0){
               #   save(D, file =  file.path(DIR_tmp, paste0('Dictionary_iter', t0+t-1, '.RData')))
               # }
             }
             t0=t0+t
             if(verbose > 1){
               save(D, file =  file.path(DIR_tmp, paste0('Dictionary_epoch', epoch, '.RData')))
               save(grad_t_1, file =  file.path(DIR_tmp, paste0('Momentum_matrix_epoch', epoch, '.RData')))
             }
           }
         },
         stop('Incorrect typeOptim. Choose: Nesterov or Momentum'))
  
  if(verbose >0){
    return(list('D' = D, 'objFunProcess' = value_OF, 
                     'learning_rate' = lr, 'Sketch' = SK_Data, 'Frequencies' = Frequencies))     
  }else{
    return(list('D' = D, 'objFunProcess' = NULL, 
                'learning_rate' = lr, 'Sketch' = SK_Data, 'Frequencies' = Frequencies))
    }
}
