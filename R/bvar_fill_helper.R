.bvar_fill_helper <- function(postdraws, tvp, nvars, tt, varname) {
  
  result <- NULL
  
  if (is.list(postdraws)) {
    
    if ("coeffs" %in% names(postdraws)) {
      if (tvp) {
        if (nvars == 1) {
          for (i in 1:tt) {
            result[[varname]][[i]] <- coda::mcmc(matrix(postdraws[["coeffs"]][(i - 1) * nvars + 1:nvars, ]))
          }
        } else {
          for (i in 1:tt) {
            result[[varname]][[i]] <- coda::mcmc(t(postdraws[["coeffs"]][(i - 1) * nvars + 1:nvars, ])) 
          } 
        }
      } else {
        result[[varname]] <- coda::mcmc(t(postdraws[["coeffs"]]))
      }        
    }
    
    if ("sigma" %in% names(postdraws)) {
      result[[paste0(varname, "_sigma")]] <- coda::mcmc(t(postdraws[["sigma"]]))
    }
    
    if ("lambda" %in% names(postdraws)) {
      result[[paste0(varname, "_lambda")]] <- coda::mcmc(t(postdraws[["lambda"]]))
    }
    
  } else {
    if (tvp) {
      if (nvars == 1) {
        for (i in 1:tt) {
          result[[varname]][[i]] <- coda::mcmc(matrix(postdraws[(i - 1) * nvars + 1:nvars, ])) 
        }
      } else {
        for (i in 1:tt) {
          result[[varname]][[i]] <- coda::mcmc(t(postdraws[(i - 1) * nvars + 1:nvars, ])) 
        } 
      }
    } else {
      result[[varname]] <- coda::mcmc(t(postdraws))     
    }
  }
  
  return(result)
}