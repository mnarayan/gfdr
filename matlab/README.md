# Description

This toolbox contains matlab implementations for a number of false discovery rate control procedures

## Usage

### Storey's Direct FDR 

The `directfdr` package contains two functions 

- `directfdr.create_options` which creates options structure
- `directfdr.run` which takes test statistics as input and returns a table of estimated fdr values and adjusted p-values 

```MATLAB
opts = directfdr.create_options()
[~,~,results] = directfdr.run(Tobs,Tperm,opts)
```

To learn more about the implementation use help

```MATLAB

help directfdr.run

  Estimate False Discovery Rate using Storey's permutation method


  Input
    t  - observed test statistic of length n_hyp x 1
    tk - permuted test statistics of length n_hyp x n_perm
    opts - (optional) defaults to create_options()

  Output
    fdr_hat     - Estimated FDR
    pval_adj    - Adjusted p-values
    results     - Tabulated results
    opts        - options structure with intermediate parameteres;
```