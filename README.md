# {boundr}: Bounded Effects and Counterfactuals for Imperfect Experiments

*This R package is under development.*

This R package is used to build and analyze structural graphical models as proposed by Pearl [2009. Causality: Models, Reasoning and Inference (Second Ed). Cambridge: Cambridge University Press. chapter 8]. This approach allows for bounded estimation of counterfactuals in an experiment where (i) all variables are discrete and (ii) not all variables are directly manipulated by the randomized assignment mechanism. This is done by modeling the correlated latent background variables determining all observable variables. Inference is conducted by constructing a Bayesian statistical model.  

To install:

```
install.packages("remotes")
remotes::install_github("karimn/boundr")
```

Please note that this package is still very incomplete, particularly in terms for documentation. 

What you can do with `boundr`:

1. Define a directed cyclic graph (DAG) composed of

        a) Discrete variables 
        b) Directed edges between variables
        c) Response functions for each variable dependent on finite types or classes (representing latent background variables)
        d) Terminal variables can be discretized continuous variables 
    
2. Define causal quantities.
3. Run full Bayesian sampling to calculate causal quantities.
