# Multi-stage Group Testing with (r,s)-regular design Algorithms

This repository simulates multi-stage (r,s)-regular design algorithms within the framework of group testing. It computes:  

- The expected, average, minimum, maximum number of tests per member  
- The expected, average, minimum, maximum duration per member  
- Additional performance criteria for evaluation such as the rate

The repository serves as a foundation for replication.  

## Technical Details  

For in-depth derivations and explanations of the group testing algorithms, refer to:  

**Balzer M. (2025).**  
*Multi-Stage Group Testing with (r,s)-Regular Design Algorithms.*  
[arXiv:2504.00611](https://arxiv.org/abs/2504.00611)  

## Example 
```
require(optimx)
require(tibble)
library(foreach)

# Set population size, probability and number of simulation
n = 1000
p = 0.05
sims = 10

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

# Computes expected, average, minimum and maximum number of tests as well as the rate of the two-stage (1,s)-regular design algorithm
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

result = dorfman(n,p,sims)

result
```

