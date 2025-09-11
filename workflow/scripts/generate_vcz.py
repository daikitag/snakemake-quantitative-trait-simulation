import numpy as np
import msprime
import pandas as pd

import tskit
import collections

import bio2zarr.tskit as ts2z

def simulate_neutral_mutation(ts, neutral_rate, seed):
    # Simulate neutral mutations on top of a SLiM tree sequence output
    # Taken from https://tskit.dev/pyslim/docs/latest/tutorial.html#adding-neutral-mutations-to-a-slim-simulation
    ts = msprime.sim_mutations(
           ts,
           rate=neutral_rate,
           keep=True,
           random_seed=seed,
    )
    return ts

def count_site_alleles(ts, tree, site):
    counts = collections.Counter({site.ancestral_state: ts.num_samples})
    for m in site.mutations:
        current_state = site.ancestral_state
        if m.parent != tskit.NULL:
            current_state = ts.mutation(m.parent).derived_state
        # Silent mutations do nothing
        if current_state != m.derived_state:
            num_samples = tree.num_samples(m.node)
            counts[m.derived_state] += num_samples
            counts[current_state] -= num_samples
    return counts

def maf_threshold(ts, maf):
    remove_site = []
    
    tree = tskit.Tree(ts)
    
    for i in range(ts.num_sites):
        site = ts.site(i)
        tree.seek(site.position)
        counts = count_site_alleles(ts, tree, site)
        # counts is a Counter object from collections
        max_allele_count = counts.most_common(1)[0][1]
        freq = max_allele_count / ts.num_samples
        if freq > (1 - maf):
            remove_site.append(i)
        # remove sites with 3+ alleles
        elif len(counts) > 2:
            remove_site.append(i)
        # remove sites with 0 frequency
        elif freq == 0:
            remove_site.append(i)

    down_sample_ts = ts.delete_sites(remove_site)
    return down_sample_ts

def main():    
    ts = tskit.load(snakemake.input.ts)

    ts = simulate_neutral_mutation(
        ts=ts, neutral_rate=float(snakemake.params.neutral_mu), seed=int(snakemake.params.neutral_seed)
    )
    
    ts = maf_threshold(ts, maf=float(snakemake.params.maf))

    individual_id_df = pd.read_csv(snakemake.input.individual_id)
   
    model_mapping = ts.map_to_vcf_model(
        individuals=individual_id_df["individual_id"],
        individual_names=individual_id_df["plink_id"]
    )

    ts2z.convert(
        ts,
        vcz_path=snakemake.output[0],
        model_mapping=model_mapping,
        worker_processes=snakemake.threads,
        contig_id=str(snakemake.params.chromosome)
    )

    
if __name__ == '__main__':
    main()