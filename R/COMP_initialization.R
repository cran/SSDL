#' @title COMP dictionary initialization 
#' @description Dictionary initialization using the Compressive Orthogonal Matching Pursuit (COMP) method  
#' @inheritParams CDL
#' @param maxIter is a maximum number of iterations in the computation of new dictionary element. 
#' The default value is 1500.
#' @param lower is a lower boundary. It is an \eqn{s}-dimensional vector.
#' @param upper is an upper boundary. It is an \eqn{s}-dimensional vector.
#' @param HardThreshold indicates whether to execute the hard thresholding step. The default is FALSE.
#' @param print_level controls how much output is shown during the optimization process. Possible values:
#' \itemize{
#'     \item 0 no output (default value)
#'     \item 1 show iteration number and value of objective function
#'     \item 2 1 + show values of weights}
#' @details The initialization routine is based on the Compressive Orthogonal Matching Pursuit 
#' (COMP) algorithm. 
#' COMP is an iterative greedy method that builds a dictionary operating on a 
#' compressed data version (a.k.a. data sketch). 
#' It alternates between expanding the dictionary \eqn{D} with a new element \eqn{d_i}, 
#' whose sketch \eqn{SK(d_i)} is the most correlated to the residue, and 
#' calculating the weights of the dictionary elements \eqn{w_1, \dots, w_K} 
#' by minimizing the difference between the data sketch \eqn{SK(Data)} and a linear combination 
#' of dictionary sketches, \emph{i.e.} \eqn{\|SK(Data) - \sum_{i=1}^K w_i \cdot SK(d_i)\|}.
#' Unlike COMP, the implemented dictionary initialization routine does not perform an additional global
#' optimization with respects to both variables: weights and dictionary elements. 
#' @return a list 
#'      \itemize{
#'      \item \code{D} is the obtained dictionary,
#'      \item \code{weights} is the resulting weights,
#'      \item \code{ObjF} is the objective function values computed at each iteration. 
#'      \item \code{Sketch} is the data sketch
#'      \item \code{Frequencies} is the frequency matrix}
#' @examples
#' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
#' lb = apply(X, 1, min)
#' ub = apply(X, 1, max)
#' X_fbm = bigstatsr::FBM(init = X, ncol = ncol(X), nrow = nrow(X))
#' m = 64
#' W = chickn::GenerateFrequencies(Data = X_fbm, m = m, N0 = ncol(X_fbm),
#'                                 ncores = 1, niter= 3, nblocks = 2, sigma_start = 0.001)$W
#' SK= chickn::Sketch(X_fbm, W)
#' D0 = COMP_initialization(K = 10, Data = X_fbm, SK_Data = SK, Frequencies = W,
#'                         lower = lb, upper = ub)$Dictionary
#' @importFrom stats optim
#' @importFrom pracma lsqnonneg lsqlincon
#' @importFrom chickn Sketch GenerateFrequencies
#' @seealso \code{\link{ObjFun_COMP_cpp}}, \code{\link{Gradient_COMP_cpp}}, 
#' \code{chickn::Sketch}, \code{chickn::GenerateFrequencies}, [chickn](https://CRAN.R-project.org/package=chickn)
#' @note COMP method has been presented in \insertRef{DBLP:journals/corr/KerivenTTG16}{SSDL} 
#' @export
COMP_initialization <- function(K,Data,SK_Data = NULL,Frequencies = NULL,
                               lower = -Inf, upper = Inf, maxIter = 1500,
                               HardThreshold = FALSE,
                               print_level = 0, ncores = 1, m = nrow(Frequencies), ...){
  
  if(length(Frequencies)==0){
    Frequencies = chickn::GenerateFrequencies(Data = Data, m = m, ...)$W
  }
  if(length(SK_Data) ==0){
    SK_Data = chickn::Sketch(Data = Data, W = Frequencies, ind.col = 1:ncol(Data),
                             ncores = ncores, parallel = FALSE)    
  }
  if(nrow(Data) != ncol(Frequencies)){
    stop('Error: incompatible data and frequency matrix dimensions')
  }
  if(length(SK_Data) != 2*nrow(Frequencies)){
    stop(paste0('Error: incompatible sketch length ', length(SK_Data), 
                ' and frequency matrix dimension', nrow(Frequencies)))
  }
  #m = nrow(Frequencies)
  d = ncol(Frequencies)
  r = SK_Data
  if(HardThreshold){
    HT_ind = 0
    iters = c(1:K, rep(K+1, K))
    D = matrix(NA, ncol = K+1, nrow = d) # K column 
    SK_D = matrix(NA, ncol = K+1, nrow = 2*m)
    prodWd = matrix(NA, ncol = m, nrow = 1)
    norm_SK_D = matrix(NA, ncol = K+1, nrow = 1)
    ObjF = rep(NA, K+1)
  }else{
    iters = 1:K
    D = matrix(NA, ncol = K, nrow = d) # K column 
    SK_D = matrix(NA, ncol = K, nrow = 2*m)
    prodWd = matrix(NA, ncol = m, nrow = 1)
    ObjF = rep(NA, K)
    
  }
  weight = rep(NA, K)
  for( t in iters){
    #=== New centroid ===
    if(print_level == 1 || print_level ==2){
      if(HardThreshold & t > K){
        cat('ITERATION', t+HT_ind, '\n')
        HT_ind = HT_ind+1
      }else{
        cat('ITERATION', t, '\n')
      }
      cat('Find new dictionary atom...')
    }
    res_D <- stats::optim(par = Data[,t], fn = ObjFun_COMP_cpp, gr = Gradient_COMP_cpp, method = 'L-BFGS-B',
                   lower = lower, upper = upper, W = Frequencies, residue = r,
                   control = list("maxit" = maxIter, 'trace' = 0))
    D[,t] = res_D$par
    #=== compute sketch of C ===
    
    prodWd = Frequencies%*%D[,t] 
    SK_D[,t] = c(cos(prodWd), sin(prodWd)) # ????Normalization /sqrt(m) 
    if(HardThreshold){
      norm_SK_D[,t] = sqrt(sum(SK_D[,t]*SK_D[,t]))
    }
    if(t > K & HardThreshold){
      if(print_level==1 || print_level ==2){cat('   Hard Thresholding...')}
      res_weight_ht = pracma::lsqnonneg(as.matrix(SK_D[,1:t], ncol = t, nrow = 2*m), as.vector(SK_Data))
      #res_weight_ht = lsqlincon(C = as.matrix(SK_D[,1:t], ncol = t, nrow = 2*m), d = as.vector(SK_Data),
      #                       A = NULL, b = NULL, Aeq = matrix(rep(1,t), ncol = t, nrow = 1), beq = as.vector(1),
      #                       lb = as.vector(rep(0, t)))
      beta = res_weight_ht$x*norm_SK_D
      beta_sort = sort.int(beta, index.return = T, decreasing = TRUE)
      
      if(max(beta_sort$ix[1:K])== K){
        if(print_level==1 || print_level == 2){cat('   NO changes\n')}
      }else{
        ind = sort(beta_sort$ix[1:K])
        D[,1:K] = D[, ind]
        SK_D[,1:K] = SK_D[,ind]
      }
      
    }
    #=== compute weights ===
    # res_weight = lsqnonneg(as.matrix(SK_D[,1:t], ncol = t, nrow = 2*m), as.vector(SK_Data))
    # weight[1:t] = res_weight$x
    if(print_level == 1 || print_level ==2){
      cat('Weight computation...\n')
    }
    res_weight = pracma::lsqnonneg(C = as.matrix(SK_D[,1:min(t, K)]),
                                   d = as.vector(SK_Data))$x
    # res_weight = pracma::lsqlincon(C = as.matrix(SK_D[,1:min(t, K)], ncol = min(t, K), nrow = 2*m),
    #                                d = as.vector(SK_Data),
    #                                A = NULL, b = NULL, Aeq = matrix(rep(1,min(t, K)), ncol = min(t, K),
    #                                                                 nrow = 1), beq = as.vector(1),
    #                                lb = as.vector(rep(0, min(t, K))))
    weight[1:min(t, K)] = res_weight
    
    #ObjF[t] = res_weight$resid.norm
    r = SK_Data - as.matrix(SK_D[,1:min(t, K)])%*%weight[1:min(t, K)]
    ObjF[t] = t(r)%*%r
    if(print_level==1 || print_level == 2){cat('Objective function: ',ObjF[t], '\n')}
    if(print_level==2){cat('Weights:',weight[1:min(t, K)], '\n')}
  }
  # sink()
  D = D[,1:K]
  return(list('Dictionary' = D, 'weight' = weight, 'ObjF' = ObjF, 
              'Sketch' = SK_Data, 'Frequencies' = Frequencies))
}