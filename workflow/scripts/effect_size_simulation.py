import numpy as np
import pandas as pd
import tskit

from snakemake.script import snakemake as snk

def obtain_mutation_df(ts, selection_scaling):
    """
    Obtain mutation dataframe from a tree sequence that is simulated in slim
    We will multiply the final selection coefficient by -2*selection_scaling, due to
    how we are modeling selection coefficient in SLiM simulation
    
    If we let `s` as the selection coefficient in the stabilizing selection model and
    model fitness of individuals in the underdominance selection model as `1`, `1 + t`,
    and `1`, we can calculate that `s = -2 * t` (see computations in the paper).
    
    In tstrait simulation, we must need the following columns in trait dataframe to
    compute genetic values:
    - site_id : Site ID from tree sequence data
    - effect_size : This will be simulated based on selection coefficient of mutations
    - causal_allele : This can be found from derived_state in mutation
    - trait_id : This will be set to 0, because we are only simulating a single trait
    """
    mutation_list = []
    for i in range(ts.num_mutations):
        mutation = ts.mutation(i)
            
        mutation_list.append(
            {
                "site_id": mutation.site,
                "selection_coeff": mutation.metadata["mutation_list"][0]["selection_coeff"],
                "causal_allele": mutation.derived_state,
                "trait_id": 0,
            }
        )
    
    mutation_df = pd.DataFrame(mutation_list)
    mutation_df["selection_coeff"] = mutation_df["selection_coeff"]*selection_scaling*(-2)

    return mutation_df

def simulate_effect(selection_coeff, n, w2, rng):
    """
    Simulates effect size from a normal distribution based on the selection coefficient
    of mutations
    n is the degree of pleiotropic effects
    """
    effect_size = rng.normal(loc = 0, scale = np.sqrt(w2/n*selection_coeff))
    return effect_size

def sim_tstrait_pleiotropy(ts, n, w, selection_scaling, seed):
    """
    ts is the tree sequence of interest
    n is the degree of pleiotropic effects
    seed is the seed that will be used in the tstrait simulation
    w2 is the stabilizing selection parameter, which will be 1 by default
    selectin_scaling is the scaling factor for the underdominance simulation model
    """
    rng = np.random.default_rng(seed=seed)
    mutation_df = obtain_mutation_df(ts, selection_scaling)
    mutation_df["effect_size"] = mutation_df.apply(lambda row: simulate_effect(row["selection_coeff"], w2=w**2, n=n, rng=rng), axis=1)
    return mutation_df


def main():    
    ts = tskit.load(snk.input[0])
    mutation_df = sim_tstrait_pleiotropy(
        ts = ts,
        n = int(snk.params.degree),
        w = float(snk.params.w),
        selection_scaling= float(snk.params.selection_scaling),
        seed = int(snk.params.seed)
    )

    mutation_df.to_csv(snk.output[0], index=False)
    

if __name__ == '__main__':
    main()
