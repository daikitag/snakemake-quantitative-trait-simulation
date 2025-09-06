import pandas as pd
import tskit
import tstrait

from snakemake.script import snakemake as snk


def simulate_phenotype(ts, mutation_df, h2, seed):
    genetic_df = tstrait.genetic_value(ts, mutation_df)
    phenotype_df = tstrait.sim_env(
        genetic_df, h2=h2, random_seed=seed
    )
    
    return phenotype_df

def main():    
    ts = tskit.load(snk.input.ts)
    mutation_df = pd.read_csv(snk.input.mutation_df)
    
    phenotype_df = simulate_phenotype(
        ts=ts, mutation_df=mutation_df, h2=float(snk.params.h2), seed=int(snk.params.seed)
    )
    
    phenotype_df.to_csv(snk.output[0], index=False)
    
if __name__ == '__main__':
    main()
