# A Robust Optimization Approach to Network Control Using Local Information Exchange


## Introduction
This repo contains all source code that reproduce the experiments in our paper 

[A Robust Optimization Approach to Network Control Using Local Information Exchange](https://arxiv.org/pdf/2405.00148) 

**We welcome any feedback and suggestions! Note that we put in maximum effort to write high quality codes. However, they may still contain bugs or not be efficient enough.**

## Prerequisites
All optimization problems are implemented in MATLAB. The implementations rely on the following third-party software: [Gurobi](https://www.gurobi.com/), [YALMIP](https://github.com/johanlofberg/YALMIP), and [RSOME](https://www.rsomerso.com/). It is necessary to install these software and add their respective directories to the MATLAB path before running the codes. We also use the [Residential Power and Battery Dataset](https://zenodo.org/records/8219786). Data preprocessing are implemented in Python. Postprocessing data is available in the repository.

## Reproducing the results
First, clone the repo

> $ git clone https://github.com/sorooshafiee/Robust_Network_Control.git

To reproduce the simulation results, you need to run m-file scripts.