import sys
import tskit
import pandas as pd

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
                "site_position": ts.site(mutation.site).position
            }
        )
    
    mutation_df = pd.DataFrame(mutation_list)
    mutation_df["selection_coeff"] = mutation_df["selection_coeff"]*selection_scaling*(-2)

    return mutation_df

def main():    
    ts = tskit.load(sys.argv[1])
    mutation_df = obtain_mutation_df(ts, selection_scaling=1e10)
    mutation_df.to_csv(sys.argv[2], index=False)
    

if __name__ == '__main__':
    main()
