# snakemake-quantitative-trait-simulation

This snakemake workflow is used to conduct the following simulations:

1. Simulating genetic trajectories by using underdominance model in SLiM
2. Simulating effect sizes by using the simulated selection coefficient in SLiM simulation
3. Computing genetic values of individuals by using the simulated effect sizes
4. Simulating environmental noise

To run this job, modify `config.yaml` in `profiles` directory to set the correct configuration for slurm to execute the code.

The simulation parameters are specified in `config.yaml` in `config` directory.