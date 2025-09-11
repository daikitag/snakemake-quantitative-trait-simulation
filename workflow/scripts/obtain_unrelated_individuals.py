import numpy as np
import pandas as pd

from snakemake.script import snakemake as snk

def obtain_unrelated_sample_df(pedigree_df):
    pedigree_sample_df = pedigree_df.sample(frac=0.25)
    parents = np.concat((pedigree_sample_df["pedigree_p1"],pedigree_sample_df["pedigree_p2"]))
    parent_counts = np.bincount(parents)
    pedigree_sample_df["parent_counts"] = pedigree_sample_df.apply(lambda row: parent_counts[row.pedigree_p1] + parent_counts[row.pedigree_p2], axis=1)

    unrelated_df = pedigree_sample_df[pedigree_sample_df["parent_counts"] == 2]
    
    return unrelated_df

def main(): 
    seed = int(snk.params.seed) + int(snk.params.chromosome)
    rng = np.random.default_rng(seed)
    pedigree_df = pd.read_csv(snk.input[0])

    unrelated_individual_df = obtain_unrelated_sample_df(pedigree_df)

    CEU_df = unrelated_individual_df[unrelated_individual_df["population"] == "CEU"]
    YRI_df = unrelated_individual_df[unrelated_individual_df["population"] == "YRI"]

    CEU_sample_df = CEU_df.sample(n=int(snk.params.ceu_number), random_state=rng)
    YRI_sample_df = YRI_df.sample(n=int(snk.params.yri_number), random_state=rng)

    final_df = pd.concat([CEU_sample_df, YRI_sample_df])
    final_df.to_csv(snk.output[0], index=False)


if __name__ == '__main__':
    main()