# Replication Package - Network-Aware Synthetic Control

This repository accompanies the study **Network-Aware Synthetic Control — Bias Correction and Regularization under Interference** von Waldow (2026), M.Sc. thesis, University of Bonn.


## Replication

1.  For replication of the studies result open `nasc_replication.Rproj` or set the working directory.
2.  Install the accompanying `nasc` R-package via `r pak::pkg_install("fvwaldow/nasc")`r and additionally `Synth`, `rstan`, `StanHeaders`, `igraph`, and `future.apply`.
3.  Run `master.R` for replication of Tables and Figures using the simulation results saved in `data/`.
4.  For a re-run of the MC simulation , run the files `simulation_BAYESIAN.R` and `simulation_SCM.R`, with `master.R` being sourced first. Both write to `output/` and never overwrite the initial results in `data/`. Even though the simulation runs in parallel, in particular the Bayesian estimations are computational heavy and take some considerable runtime.


## Folder Structure

```         
.
├── master.R                   
├── README.md
├── nasc_replication.Rproj
├── R/
│   ├── dgp_functions.R                   data-generating process
│   ├── MC Simulations/
│   │   ├── simulation_BAYESIAN.R         MC: BSCM / BC / CR / NASC
│   │   └── simulation_SCM.R              MC: SCM
│   ├── Figures/
│   │   ├── Plot_Ternary_ETDir.R          ETDir simplex geometry
│   │   ├── Plot_MC_bias.R                bias boxplots
│   │   ├── Plot_MC_SC_contam.R           estimated contamination boxplots
│   │   ├── Plot_MC_bias_corr_factor.R    bias-correction factor boxplots
│   │   └── Plot_example_MC_iter.R        single-replicate illustration
│   └── Tables/
│       ├── Table_MC_results.R            bias / RMSE / coverage
│       └── Table_MC_SC_contam.R          estimated contamination
├── data/
│   ├── mc_result_BAYESIAN.rds            MC: BSCM / BC / CR / NASC
│   ├── mc_result_SCM.rds                 MC: SCM
│   └── networks/                         adjacency matrices
└── output/                     
    ├── figures/
    ├── tables/
    └── sessionInfo.txt                   software environment
```



All R files were run on R version 4.4.2 on Windows 11 x64. 

Frederik von Waldow
