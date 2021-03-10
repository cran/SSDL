#' Additional function for the TV norm initialization procedure.
#' @inheritParams make_copy
#' @param step is a step size. The default value is 1.
#' @return three dimensional vector with the optimal step size, lower and upper step size bounds.
#' @keywords internal
#' @export 
step <- function(x, n_copies, step = 1){
  while(step < x['rb']){
    n_current = floor((x['rb'] - x['x2'])/step) + floor((x['x1'] - x['lb'])/step)
    if(n_current == n_copies -1){
      break
    }else if(n_current > n_copies - 1){
      step = step + 1
    }else{
      step = step -1
      break
    }
  }
  if(step == 0){
    stop(paste('impossible make',n_copies, 'signal copies', sep = ' '))
  }else{
    return(c(step, floor((x['x1'] - x['lb'])/step), floor((x['rb'] - x['x2'])/step)))
  }
}

#' function makes copy of the selected signal
#' @param x is a selected signal pattern
#' @param n_copies is a number of the copies.
#' @return a matrix, which contains \code{n_copies} of the patters signal in the columns.
#' @keywords internal
#' @export
make_copy <- function(x, n_copies = 2){
  step = x[length(x)-4]
  kl = x[length(x)-3]
  kr = x[length(x)-2]
  x1 = x[length(x)-1]
  x2 = x[length(x)]
  if(kl + kr > n_copies){
    if(kl < kr){
      d = round((kl+kr)/2) - kl
      kl = round(n_copies/2)-1 - d 
      kr = round(n_copies/2) + d
    }else{
      d = round((kl+kr)/2) - kr
      kr = round(n_copies/2)-1 - d
      kl = round(n_copies/2) + d
    }
    
  }
  data = matrix(NA, nrow = (length(x)-5), ncol = kl+kr)
  
  k = 1
  while(k < kl+1){
   
    data[,k] = c(x[(1+step*k+1): (length(x)-5)], x[1:(1 + step*k)])
    k = k+1
  }
  #data[,k] = x[1:(length(x)-5)]
  k = 1
  while(k < kr + 1){
    #data[,k + kl] = c(x[(x2 + step*k): (length(x)-5)], x[1:(x2 + step*k -1)])
    data[,k + kl] = c(x[(length(x)-5- step*k+1): (length(x)-5)], x[1:(length(x)-5- step*k)])
    k = k+1
  }
  return(data)
}
#' @title TV norm dictionary initialization
#' @description Dictionary initialization using a TV norm criterion
#' @inheritParams CDL
#' @param cutoff is a cut off value, the default value is 0.5.
#' @param Npattern is a number of patterns selected in the dataset to create the dictionary 
#' @param set_size is a maximum size of the set of possible patterns. 
#' @param DoCopies indicates whether to duplicate patterns.
# @param type is a type of the criterion of the pattern selection. It can be either the correlation ('Corr), 
# the random selection ('rand) or Wasserstein distance ('W_dist).
#' @param ncores is a number of cores
#' @param DIR_tmp is a directory to save temporary files 
#' @return a dictionary matrix 
#' @importFrom bigstatsr big_copy big_apply
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @importFrom utils tail head
#' @details The dictionary is initialized by extracting and duplicating patterns with the highest TV norm values
#' To limit the set of possible patterns, only signals with the correlation less then a fixed threshold \code{cutoff} are 
#' taken into account. If the set of possible patterns is too large, it can be further reduced by taking only 
#' \code{set_size} less correlated patterns. The implemented initialization routine can only be applied to positive value data.   
#' @examples
#' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
#' X_fbm = bigstatsr::FBM(init = X, ncol = ncol(X), nrow = nrow(X))
#' D0 = TV_initialization(X_fbm, K = 20, Npattern = 5, DoCopies = TRUE, ncores = 1)
#' @export
TV_initialization<- function(Data, K, cutoff = 0.5, Npattern = 8,  set_size = ncol(Data),  
                        DoCopies = FALSE, 
                        #type = 'Corr', 
                        ncores = 4,
                        DIR_tmp = tempdir()){
  # Normalize data =======================================================
  bf_tvmax = paste0(DIR_tmp, 'Data_tvmax')
  if(file.exists(paste0(bf_tvmax, '.bk'))){
    file.remove(paste0(bf_tvmax, '.bk'))
  }
  if(file.exists(paste0(bf_tvmax, '.rds'))){
    file.remove(paste0(bf_tvmax, '.rds'))
  }
  # make a copy of data 
  Data_tvmax <- bigstatsr::big_copy(Data, backingfile = bf_tvmax)$save()
  # normalize by max value
  bigstatsr::big_apply(Data_tvmax, a.FUN = function(X, ind){
    X[, ind] <- X[, ind]/matrix(apply(X[,ind], 2, max), nrow = nrow(X), ncol = ncol(X[,ind]), 
                                byrow = TRUE)
    NULL
  }, a.combine = 'cbind', ncores = ncores)
  # compute TV norm in the matrix columns
  tv <- bigstatsr::big_apply(Data_tvmax, a.FUN = function(X, ind){
    colSums(abs(X[2:nrow(X),ind] - X[1:(nrow(X)- 1), ind]))
  }, a.combine = 'c', ncores = ncores)
  # tv <- big_apply(Data_tvmax, a.FUN = function(X, ind){
  #   apply(X[,ind], 2, max)
  # }, a.combine = 'c',ncores = ncores)
  
  # Take not correlated signals with maximum TV norm  =====
  ind_tvmax = sort.int(tv, decreasing = FALSE, index.return = TRUE)$ix[1:set_size]
  if(file.exists(paste0(bf_tvmax, 'notcorr.bk'))){
    file.remove(paste0(bf_tvmax, 'notcorr.bk'))
  }
  if(file.exists(paste0(bf_tvmax, 'notcorr.rds'))){
    file.remove(paste0(bf_tvmax, 'notcorr.rds'))
  }
  Data_tvmax_notcorr <- bigstatsr::big_copy(Data_tvmax, ind.col = ind_tvmax, 
                                  backingfile = paste0(bf_tvmax, 'notcorr'))$save()
  df_tvmax = as.data.frame(Data_tvmax_notcorr[])
  vars = paste0('Signal', 1:ncol(Data_tvmax_notcorr[]))
  names(df_tvmax) = vars
  
  mcor = 1-abs(cor(df_tvmax))
  df <- reshape2::melt(mcor)
  
  # switch(type,
  #        'Corr'= {
  #          mcor = 1-abs(cor(df_tvmax))
  #          df <- reshape2::melt(mcor)
  #        },
  #        'W_dist'={
  #          cumsum_tvmax = apply(df_tvmax, 2, cumsum)
  #          last = cumsum_tvmax[nrow(cumsum_tvmax),]
  #          cumsum_tvmax = cumsum_tvmax/ last[col(cumsum_tvmax)]
  #          W_dist = matrix(NA, ncol = ncol(Data_tvmax_notcorr), nrow = ncol(Data_tvmax_notcorr))
  #          for(col_current in 1:ncol(Data_tvmax_notcorr)){
  #            for(j in 1:ncol(Data_tvmax_notcorr)){
  #              W_dist[col_current,j] = sum(abs((cumsum_tvmax[,col_current] - cumsum_tvmax[,j])))
  #            }
  #          }
  #          colnames(W_dist) = vars
  #          rownames(W_dist) = vars
  #          df <- reshape2::melt(W_dist)
  #        },
  #        'rand' = {
  #          ind.dict = sample(1:ncol(df_tvmax), K)
  #          df <- df_tvmax[,ind.dict]
  #        },
  #        stop("Unkown type value"))

    if(!DoCopies){
      D = as.matrix(df)
    }else{
      if(K%%Npattern!= 0){
        n_copies = K%/%Npattern
      }else{
        n_copies = K / Npattern
      }
      if(n_copies > nrow(Data_tvmax_notcorr)){
        warning('number of copies bigger than number of discrete points')
      }
   
      names(df) = c('Signal_row', 'Signal_col', 'value')
      Signal_row<-Signal_col<-value <- NULL
      # Choose Npattern for duplicate
      df_cutoff = df[which(df$value > cutoff),]
      unique_vars = unique(df_cutoff$Signal_col)
      pattern = vector(mode = 'character')
      #pattern[1] = as.character(unique(df_cutoff$Var2))[1]
      pattern = as.character(unique(df_cutoff$Signal_col))[1:Npattern]
      k=1
      set_candidate = vars
      best_pattern = pattern
     
      while(length(pattern) < Npattern & k <= length(unique_vars)){
        
        df_current = subset(df_cutoff, Signal_col == tail(pattern, 1) & !(Signal_row %in% pattern)
                            & (Signal_row %in% set_candidate))
        
        if(length(intersect(set_candidate, df_current$Signal_row))!=0){
          new = as.character(df_current$Signal_row[which.max(abs(df_current$value))])
          pattern = c(pattern,new)
          set_candidate = setdiff(intersect(set_candidate, df_current$Signal_row), pattern)
        }else{
          k=k+1
          if(k <= length(unique_vars)){
            pattern = as.character(unique_vars)[k]
          }
          set_candidate = vars
        }
        if(length(pattern) > length(best_pattern)){
          best_pattern = pattern
        }
      }
      df_pattern = as.data.frame(Data_tvmax[, which(vars %in% best_pattern)])
      names(df_pattern) = best_pattern
    #  df_pattern = as.data.frame(apply(df_pattern, 2, function(x) x/sqrt(sum(x*x))))
      D = matrix(NA, nrow= nrow(Data_tvmax), ncol = K)
      if(n_copies ==1){
        for(i in 1:ncol(df_pattern)){
          D[,i]= df_pattern[,i]
        }
       
        if(ncol(df_pattern) < K){
          ind_patterns = sample(ncol(df_pattern), K - ncol(df_pattern))
          D[, (ncol(df_pattern)+1):K] = as.matrix(df_pattern[,ind_patterns])
          
        }
      }else{
        i=1
        for(x in df_pattern){
          D[,i] = x
          i = i+1
          cs = cumsum(x)
          cs = cs/cs[length(cs)]
          interval = c(max(min(which(cs > 0))-2,1), min(max(which(cs < 1))+1, 
                                                        nrow(Data_tvmax_notcorr)))
          signal_width = interval[2] - interval[1]
          
          Grid_size = nrow(Data_tvmax_notcorr)%/%(n_copies-1)
          diff = nrow(Data_tvmax_notcorr) - (n_copies-1)*Grid_size - signal_width
          signal_copy = c(tail(x, -interval[1]), head(x, interval[1]))
          D[,i] = signal_copy
          i=i+1
          k = 1
          
          while(k < n_copies - 1){
            if(k > n_copies -1 + diff & diff < 0 ){
              if(Grid_size == 1){
                break()
              }else{
                Grid_size = Grid_size - 1 
                diff = 0
              }
              
            }
            signal_copy = c(tail(signal_copy, Grid_size), head(signal_copy, - Grid_size))
            D[,i] = signal_copy
            i=i+1
            k = k+1
          }
        }  
        if(i < K){
          ind_patterns = sample(ncol(df_pattern), K - i+1)
          D[, i:K] = as.matrix(df_pattern[,ind_patterns])
        }
      }
    }
  if(file.exists(paste0(bf_tvmax, '.bk'))){
    file.remove(paste0(bf_tvmax, '.bk'))
  }
  if(file.exists(paste0(bf_tvmax, '.rds'))){
    file.remove(paste0(bf_tvmax, '.rds'))
  }
  if(file.exists(paste0(bf_tvmax, 'notcorr.bk'))){
    file.remove(paste0(bf_tvmax, 'notcorr.bk'))
  }
  if(file.exists(paste0(bf_tvmax, 'notcorr.rds'))){
    file.remove(paste0(bf_tvmax, 'notcorr.rds'))
  }
  return(D)
}
