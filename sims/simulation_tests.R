library(optimx)
library(tidyverse)
library(doParallel)
library(foreach)
library(vroom)

# Parallel backend
registerDoParallel(makeCluster(detectCores()))

calculate = function(n, p, sims = 0) {
  
# Binary entropy function
bound = function(p) {
  if(round(n*p) != 0) {
    res = -p * log2(p) - (1-p)*log2(1-p)
    return(res)
  } else {
    res = NaN
    return(res)
  }
}
  
  
# Two-Stage
dorfman = function(n, p, sims = 0) {
  
  # Searching for the No. Tests
  if(round(n * p) == 0) {
    theo = 1
    opts = 1
    num_tests = 1
  } else {
    opt = function(p, s) {
      res = ((1/s) + 1 - (1 - p)^s)
    }
    
    optimization = optimx(par = c(s = 1), fn = function(params) opt(p, params["s"]), method = "L-BFGS-B")
    
    if(optimization$s > n | optimization$value > 1 | optimization$convcode == 1) {
      theo = n
      opts = 1
    } else {
      theo = optimization$value * n
      opts = optimization$s
    }
  }
  
  # If sims != 0 simulate the procedure
  if (sims != 0) {
    if (round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      # Use foreach to parallelize the loop
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        set.seed(s)
        # Simulate the Dorfman procedure
        if ((opts - 1) != 0) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          num_groups = ceiling(n / opts)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts))
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          
          
          # Stage 2: Test individuals in positive groups individually
          num_tests = num_groups + length(unlist(groups[p_groups], use.names = FALSE))
          
        } else {
          num_tests = n
        }
        
        # Return num_tests value for this iteration
        return(num_tests)
      }
      
      
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA
  }
  
  
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) = "SP-Two"
  
  return(df)
}

grid = function(n, p, sims = 0) {
  
  # Searching for the No. Tests
  if(round(n * p) == 0) {
    theo = 1
    opts = 1
    num_tests = 1
  } else {
      opt = function(p, s) {
        res =  ((2/s) + p + (1 - p) * (1 - (1 - p)^(s - 1))^2)
      }
      
      optimization = optimx(par = c(s = 1), fn = function(params) opt(p, params["s"]), method = "L-BFGS-B")
      
      if(optimization$s > n | optimization$value > 1 | optimization$convcode == 1) {
        theo = n
        opts = 1
      } else {
        theo = optimization$value * n
        opts = optimization$s
      }
    }
  
  # If sims != 0 simulate the procedure
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        set.seed(s)
        if ((opts - 1) != 0) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          shuffled_indices2 = sample(n)
          num_groups = ceiling(n / opts)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts))
          groups2 = split(shuffled_indices2, ceiling(seq_along(1:n)/opts))
          
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          p_groups2 = which(sapply(groups2, function(g) any(g %in% infected)))
          
          mts = intersect(unlist(groups[p_groups], use.names = FALSE), unlist(groups2[p_groups2], use.names = FALSE))
          
          
          # Stage 2: Test individuals in positive groups individually
          num_tests = 2*num_groups + length(mts)
          
        } else {
          num_tests = n
        }
        return(num_tests)
      }
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA
  }
  


  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) = "DP-Two"
  
  return(df)
}

rpooling = function (n, p, sims = 0) {
  
  # Optimization of rpooling
  
  if(round(n * p) == 0) {
    theo = 1
    opts = 1
    num_tests = 1
  } else {
    
    opt = function(p, s, r) {
      res =  ((r/s) + p + (1 - p) * (1 - (1 - p)^(s - 1))^r)
    }
    
    optimization = optimx(par = c(s = 1, r = 1), fn = function(params) opt(p, params["s"], params["r"]), lower = c(1,1), method = "L-BFGS-B")
    
    if(optimization$s > 1000 | optimization$r < 1 | optimization$value > 1 | optimization$convcode == 1) {
      theo = n
      opts = 0
      optr = 0
    } else {
      theo = optimization$value * n
      opts = optimization$s
      optr = round(optimization$r)
    }
  }
  
  # If sims != 0 simulate the procedure
  
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = numeric(sims)  # Initialize a vector to store num_tests values
      
      for (s in 1:sims) {
        if ((opts - 1) > 0 & optr >= 1) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          
          groups_list = list()
          
          # GeneRate and split shuffled indices for each set
          num_groups = ceiling(n / opts)
          
          for (i in 1:optr) {
            shuffled_indices = sample(n)
            groups = split(shuffled_indices, ceiling(seq_along(1:n) / opts))
            groups_list[[i]] = groups
          }
          
          
          p_groups_list = list()  # Initialize vector to store indices of positive groups
          
          
          for(i in 1:optr) {
            p_groups = c()
            for (j in 1:length(groups_list[[i]])) {
              if (sum(groups_list[[i]][[j]] %in% infected) > 0) {
                # If group has at least one infected individual, save its index
                p_groups = c(p_groups, j)
              }
              p_groups_list[[i]] = p_groups
            }
          }
          
          
          # Stage 2: Test individuals in positive groups individually
          ps = list()
          for(i in 1:optr) {
            ps[[i]] = unlist(unname(groups_list[[i]][p_groups_list[[i]]]), use.names = FALSE)
          }
          
        
          # IteRate through each list starting from the second one
          
          if(optr == 1) {
            common_elements = ps[[1]]
          } else {
            common_elements = ps[[1]]
            for (i in 2:length(ps)) {
              # Find the intersection between the current common elements and the next list
              common_elements = intersect(common_elements, ps[[i]])
            }
          }
          
          num_tests = optr*num_groups + length(common_elements)
          
          
        } else {
          num_tests = n
        }
        # Save num_tests value for this iteration
        num_tests_vector[s] = num_tests
      }   
      
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA
  }
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) = "RP-Two"
  
  return(df)
}

# Three-Stage
three = function(n, p, sims = 0) {
  
  # Searching for the No. Tests
  if(round(n * p) == 0) {
    theo = 1
    opts1 = 0
    opts2 = 0
  } else {
    opt =  function(p, s1, s2) {
      res = (1/s1 + 1/s2*(1 - (1-p)^s1) + (1 - (1-p)^s2))
    }
    
    optimization = optimx(par = c(s1 = 1, s2 = 1), fn = function(params) opt(p, params["s1"], params["s2"] ), method = c("L-BFGS-B"), lower = c(1,1))
    
    
    if(optimization$s1 > n | optimization$value > n |  optimization$kkt2 == FALSE | optimization$value == -Inf | optimization$convcode == 1) {
      theo = n
      opts1 = 0
      opts2 = 0
    } else {
      theo = optimization$value * n
      opts1 = optimization$s1
      opts2 = optimization$s2
    }
  }
  
  
  # Simulate Procedure 
  
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        set.seed(s)   
        if (opts1 != 0 & opts2 != 0) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          num_groups = ceiling(n / opts1)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts1))
          
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          
          # Stage 2: Divide positive pools in subpools
          n2 = as.numeric(length(unlist(groups[p_groups], use.names = FALSE)))
          shuffled_indices2 = unname((unlist(groups[p_groups], use.names = FALSE)))
          num_groups2 = ceiling(n2/opts2)
          groups2 = split(shuffled_indices2, ceiling(seq_along(1:n2)/opts2))
          
          
          p_groups2 = which(sapply(groups2, function(g) any(g %in% infected)))
          
          # Stage 3: Test individuals in positive groups individually
          num_tests = num_groups + num_groups2 + length(unlist(groups2[p_groups2], use.names = FALSE))
          
      
          
        } else {
          num_tests = n
        }
        
        # Return both num_tests and duration for this iteration
       return(num_tests)
      }
      # Calculate statistics
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)

      
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA

  }
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) = "SP-Three"
  
  return(df)
}

onethree = function(n, p, sims = 0) {
  
  if(round(n * p) == 0) {
    theo = 1
    opts1 = 0
    opts2 = 0
  } else {
    opt =  function(p, s1, s2) {
      res = (2/s1 + 1/s2 * (p + (1-p)* (1 - (1 - p)^(s1 - 1))^2) + (1 - (1-p)^s2)) 
    }
    
    optimization = optimx(par = c(s1 = 1, s2 = 1), fn = function(params) opt(p, params["s1"], params["s2"] ), lower = c(1,1), method = c("L-BFGS-B"))
    
    
    if(optimization$s1 > n | optimization$value > n | optimization$kkt1 == FALSE | optimization$value == -Inf | optimization$convcode == 1) {
      theo = n
      opts1 = 0
      opts2 = 0
    } else {
      theo = optimization$value*n
      opts1 = optimization$s1
      opts2 = optimization$s2
    }
  }
  
  
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        if (opts1 != 0 & opts2 != 0) {
          set.seed(s)
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          shuffled_indices2 = sample(n)
          num_groups = ceiling(n / opts1)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts1))
          groups2 = split(shuffled_indices2, ceiling(seq_along(1:n)/opts1))
          
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          p_groups2 = which(sapply(groups2, function(g) any(g %in% infected)))
          
          mts = intersect(unlist(groups[p_groups], use.names = FALSE), unlist(groups2[p_groups2], use.names = FALSE))
          
          
          
          # Stage 2: Divide positive pools in subpools
          n2 = as.numeric(length(mts))
          shuffled_indices3 = mts
          num_groups2 = ceiling(n2/opts2)
          groups3 = split(shuffled_indices3, ceiling(seq_along(1:n2)/opts2))
          
          
          p_groups3 = which(sapply(groups3, function(g) any(g %in% infected)))
          
  
          
          
          # Stage 3: Test individuals in positive groups individually
          num_tests = 2*num_groups + num_groups2 + length(unlist(groups3[p_groups3], use.names = FALSE))
          
          
          
        } else {
          num_tests = n
        }
        
        # Return both num_tests and duration for this iteration
        return(num_tests)
      }
      # Calculate statistics
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
      
      
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA
    
  }
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) = "DP-Three"
  
  return(df)
  
  
}

# Four-Stage
four = function(n, p, sims = 0) {
  
  if(round(n * p) == 0) {
    theo = 1
    opts1 = 0
    opts2 = 0
    opts3 = 0
  } else {
    opt =  function(p, s1, s2, s3) {
      res = (1/s1 + 1/s2*(1 - (1-p)^s1) + 1/s3*(1 - (1-p)^s2) + (1-(1-p)^s3))
    }
    
    optimization = optimx(par = c(s1 = 1, s2 = 1, s3 = 1), fn = function(params) opt(p, params["s1"], params["s2"], params["s3"]), method = c("L-BFGS-B"), lower = c(1,1,1))
    
    if(optimization$s1 > 1000 | optimization$value > n |  optimization$kkt2 == FALSE | optimization$value == -Inf | optimization$convcode == 1) {
      theo = n
      opts1 = 0
      opts2 = 0
      opts3 = 0
    } else {
      theo = optimization$value * n
      opts1 = optimization$s1
      opts2 = optimization$s2
      opts3 = optimization$s3
    }
  }
  
  
  
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        set.seed(s) 
          
        if (opts1 > 0 & opts2 > 0 & opts3 > 0) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          num_groups = ceiling(n / opts1)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts1))
          
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          
          # Stage 2: Divide positive pools in subpools
          n2 = as.numeric(length(unlist(groups[p_groups], use.names = FALSE)))
          shuffled_indices2 = unname(unlist(groups[p_groups], use.names = FALSE))
          num_groups2 = ceiling(n2/opts2)
          groups2 = split(shuffled_indices2, ceiling(seq_along(1:n2)/opts2))
          
          
          p_groups2 = which(sapply(groups2, function(g) any(g %in% infected)))
          
          # Stage 3: Divide positive pools in subpools
          n3 = as.numeric(length(unlist(groups2[p_groups2], use.names = FALSE)))
          shuffled_indices3 = unname(unlist(groups2[p_groups2], use.names = FALSE))
          num_groups3 = ceiling(n3/opts3)
          groups3 = split(shuffled_indices3, ceiling(seq_along(1:n3)/opts3))
          
          
          p_groups3 = which(sapply(groups3, function(g) any(g %in% infected)))
          
          
          # Stage 4: Test individuals in positive groups individually
          num_tests = num_groups + num_groups2 + num_groups3 + length(unlist(groups3[p_groups3], use.names = FALSE))
          
          
        } else {
          num_tests = n
        }
        
        # Return both num_tests and duration for this iteration
        return(num_tests)
      }
      
      
      # Calculate statistics
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
      
    
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA

  }
  
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) =  "SP-Four"
  
  return(df)
}

onefour = function(n, p, sims = 0) {
  
  if(round(n * p) == 0) {
    theo = 1
    opts1 = 0
    opts2 = 0
    opts3 = 0
  } else {
    opt =  function(p, s1, s2, s3) {
      res = (2/s1 + 1/s2*(p + (1-p)*(1 - (1-p)^(s1-1))^2) + 1/s3*(1 - (1-p)^s2) + (1 - (1-p)^s3))
      return(res)
    }
    
    optimization = optimx(par = c(s1 = 1, s2 = 1, s3 = 1), fn = function(params) opt(p, params["s1"], params["s2"], params["s3"]), lower = c(1,1,1), method = c("L-BFGS-B"))
    
    if(optimization$s1 > 1000 | optimization$value > 1 |  optimization$kkt2 == FALSE | optimization$value == -Inf | optimization$convcode == 1) {
      theo = n
      opts1 = 0
      opts2 = 0
      opts3 = 0
    } else {
      theo = optimization$value * n
      opts1 = optimization$s1
      opts2 = optimization$s2
      opts3 = optimization$s3
    }
  }
  
  
  if (sims != 0) {
    if(round(n * p) == 0) {
      mtests = 1
      ltests = 1
      utests = 1
    } else {
      
      num_tests_vector = foreach(s = 1:sims, .combine = c) %dopar% {
        set.seed(s) 
        
        if (opts1 > 0 & opts2 > 0 & opts3 > 0) {
          # State infected individuals
          infected = sample(n, size = round(p * n))
          
          
          # Stage 1: Divide population into random groups of size s
          shuffled_indices = sample(n)
          shuffled_indices2 = sample(n)
          num_groups = ceiling(n / opts1)
          groups = split(shuffled_indices, ceiling(seq_along(1:n)/opts1))
          groups2 = split(shuffled_indices2, ceiling(seq_along(1:n)/opts1))
          
          
          p_groups = which(sapply(groups, function(g) any(g %in% infected)))
          p_groups2 = which(sapply(groups2, function(g) any(g %in% infected)))
          
          mts = intersect(unlist(groups[p_groups], use.names = FALSE), unlist(groups2[p_groups2], use.names = FALSE))
          
          
          # Stage 2: Divide positive pools in subpools
          n2 = length(mts)
          shuffled_indices3 = mts
          num_groups2 = ceiling(n2/opts2)
          groups3 = split(shuffled_indices3, ceiling(seq_along(1:n2)/opts2))
          
          
          p_groups3 = which(sapply(groups3, function(g) any(g %in% infected)))
          
          
          # Stage 3: Divide positive pools in subpools
          n3 = as.numeric(length(unlist(groups3[p_groups3], use.names = FALSE)))
          shuffled_indices4 = unname(unlist(groups3[p_groups3], use.names = FALSE))
          num_groups3 = ceiling(n3/opts3)
          groups4 = split(shuffled_indices4, ceiling(seq_along(1:n3)/opts3))
          
          
          p_groups4 = which(sapply(groups4, function(g) any(g %in% infected)))
          
  
          
          # Stage 4: Test individuals in positive groups individually
          num_tests = 2*num_groups + num_groups2 + num_groups3 + length(unlist(groups4[p_groups4], use.names = FALSE))
       
        
        } else {
          num_tests = n
        }
        
        # Return both num_tests and duration for this iteration
        return(num_tests)
      }
      
      
      # Calculate statistics
      mtests = mean(num_tests_vector)
      ltests = min(num_tests_vector)
      utests = max(num_tests_vector)
      
      
    }
  } else {
    mtests = NA
    ltests = NA
    utests = NA
    
  }
  
  
  
  df = data.frame("n" = n,
                  "p" = p,
                  "Theoretical" = theo / n,
                  "Tests" = mtests / n,
                  "Lower" = ltests / n,
                  "Upper" = utests / n,
                  "Rate" = n * bound(p) / theo)
  
  row.names(df) =  "DP-Four"
  
  return(df)
}

Tests = rbind(dorfman(n, p, sims),
              grid(n, p, sims),
              rpooling(n, p, sims),
              three(n, p, sims),
              onethree(n, p, sims),
              four(n, p, sims),
              onefour(n, p, sims))

Tests = tibble::rownames_to_column(Tests, "Algorithm")

return(Tests)

}


n = 1000
p_val = seq(0,0.35,0.001)


results = lapply(p_val, function(p) calculate(n, p, sims = 100))

# Combine list of data frames into a single data frame
combined_results = do.call(rbind, results)

# Create a new column for the p values
combined_results$p = rep(p_val, each = nrow(results[[1]]))

#combined_results = vroom("D:/Universität/PhD/Project 2/100.csv")
filtered_results = combined_results %>% filter(Algorithm != "RP-Two")
filtered_results = filtered_results %>% filter(Theoretical != 1)

combined_results$Algorithm = factor(combined_results$Algorithm, levels = c("SP-Two", "DP-Two", "RP-Two", "SP-Three", "DP-Three", "SP-Four", "DP-Four"))
filtered_results$Algorithm = factor(filtered_results$Algorithm, levels = c("SP-Two", "DP-Two", "RP-Two", "SP-Three", "DP-Three", "SP-Four", "DP-Four"))

# filter results for evaluation of performance
pereval = combined_results %>% filter(Theoretical != 1)
pereval = pereval %>% filter(p < 0.077)

# Calculate the absolute difference in simulations study
filtered_results$Diff = 100*abs((filtered_results$Theoretical  - filtered_results$Tests)/ filtered_results$Theoretical)

res = filtered_results %>% filter(p < 0.077)
res = filtered_results %>% filter(p > 0.077 & p < 0.182)
res = filtered_results %>% filter(p > 0.182)

error = res %>%
  group_by(Algorithm) %>%
  summarise(MAPE = mean(Diff, na.rm = TRUE))

error


x11()
# Plotting results for the simulation study
ggplot(filtered_results, aes(x = p)) +
  geom_line(aes(y = Tests, linetype = "Average"), linewidth = 1, alpha = 0.75) +
  geom_line(aes(y = Theoretical, linetype = "Expectation"), linewidth = 1, alpha = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "Range"), alpha = 0.5) +
  facet_wrap(~ Algorithm, scales = "free") +
  scale_linetype_manual(name = "Legend", values = c("Average" = "solid", "Expectation" = "dashed")) +
  scale_fill_manual(name = "Uncertainty", values = c("Range" = "grey80")) +
  labs(x = "Probability", y = "Tests per member") +
  theme_bw() +
  theme(legend.position = "bottom")




# Plotting the results for the performance evaluation

ggplot(pereval, aes(x = p, y = Theoretical, color = Algorithm, linetype = Algorithm)) +
  geom_line(linewidth = 1.0) +
  scale_linetype_manual(values = c("SP-Two" = "solid",
                                   "DP-Two" = "solid",
                                   "RP-Two" = "dashed",
                                   "SP-Three" = "solid",
                                   "DP-Three" = "solid",
                                   "SP-Four" = "solid",
                                   "DP-Four" = "solid")) +
  scale_color_manual(values = c("SP-Two" = "blue",
                                "DP-Two" = "cyan",
                                "RP-Two" = "red",
                                "SP-Three" = "green",
                                "DP-Three" = "chartreuse4",
                                "SP-Four" = "darkgoldenrod1",
                                "DP-Four" = "darkgoldenrod")) +
  labs(x = "Probability", y = "Expected number of tests per member") +
  theme_bw() +
  theme(legend.title = element_blank(), 
         legend.position = "bottom",
         legend.key.width = unit(1.5, "cm")) 

# ggplot(pereval, aes(x = p, y = Theoretical, color = Algorithm, linetype = Algorithm)) + 
#   geom_line(linewidth = 1.0) +
#   scale_color_manual(values = c("SP-Two" = "#808080",   
#                                 "DP-Two" = "#4D4D4D",   
#                                 "RP-Two" = "#000000",   
#                                 "SP-Three" = "#B3B3B3", 
#                                 "DP-Three" = "#D9D9D9", 
#                                 "SP-Four" = "#666666",  
#                                 "DP-Four" = "#999999")) +
#   scale_linetype_manual(values = c("SP-Two" = "dashed",      
#                                    "DP-Two" = "twodash",     
#                                    "RP-Two" = "solid",     
#                                    "SP-Three" = "12345678",  
#                                    "DP-Three" = "dotted", 
#                                    "SP-Four" = "longdash",  
#                                    "DP-Four" = "dotdash")) +
#   labs(x = "Probability", y = "Expected number of tests per member", 
#        title = "Expected number of tests per member by algorithm") +
#   theme_bw() +
#   theme(legend.title = element_blank(), 
#         legend.position = "bottom",
#         legend.key.width = unit(1.5, "cm")) 
# 
# 
# 
ggplot(pereval, aes(x = p, y = Rate, color = Algorithm, linetype = Algorithm)) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c("SP-Two" = "solid",
                                   "DP-Two" = "solid",
                                   "RP-Two" = "dashed",
                                   "SP-Three" = "solid",
                                   "DP-Three" = "solid",
                                   "SP-Four" = "solid",
                                   "DP-Four" = "solid")) +
  scale_color_manual(values = c("SP-Two" = "blue",
                                "DP-Two" = "cyan",
                                "RP-Two" = "red",
                                "SP-Three" = "green",
                                "DP-Three" = "chartreuse4",
                                "SP-Four" = "darkgoldenrod1",
                                "DP-Four" = "darkgoldenrod")) +
  labs(x = "Probability", y = "Rate", title = "Rate for different probabilities by algorithm") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"))






