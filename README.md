# UOT

Reproducible code for the numerical experiments in the paper:

**“On statistical properties of matching via multimarginal unbalanced optimal transport”**  
Florian Gunsilius and Yuliang Xu

## Repository structure

- `./Lalonde/`  
  Code to reproduce case study results on the Lalonde dataset (ET submission version).
  
- `./simulation/`  
  Code to reproduce simulation results in Section 4 (ET submission version).

## Main scripts

- `analyze_allATO.R`  
  Runs the simulation code for the earlier arXiv version of the paper (Case 1: 2D covariates + linear outcome):  
  https://arxiv.org/pdf/2112.04398

  You can modify the simulation setup (case number, covariate dimension, outcome function) by editing:
  - `ATO.R`
