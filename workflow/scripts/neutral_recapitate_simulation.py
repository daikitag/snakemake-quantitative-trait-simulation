import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit

import collections

from snakemake.script import snakemake as snk

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

def recapitate(ts, ancestral_Ne, ratemap, seed):
    # This function conducts recapitation based on a single population with Ne
    # being ancestral_Ne. ratemap is used for recombination rate.
    ts = pyslim.recapitate(
        ts, ancestral_Ne=ancestral_Ne, recombination_rate=ratemap, random_seed=seed
        )
    return ts

def obtain_ratemap(filename):
    # Obtain genetic map from a text file with tab delimiter, and convert it to
    # msprime RateMap object
    genetic_map = pd.read_csv(filename, delimiter="\t")
    position = np.array(genetic_map["Position(bp)"])
    position = np.insert(position, 0, 0)
    
    # Rates are in cM/Mb, so we need to convert it to recombination rate
    # by multiplying 1e-8 to it
    rate = np.array(genetic_map["Rate(cM/Mb)"]) * 1e-8
    
    return msprime.RateMap(position=position, rate=rate)

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
    
    ts = tskit.load(snk.input[0])
    ratemap = obtain_ratemap(snk.params.genetic_map)
    
    ts = recapitate(
        ts=ts, ancestral_Ne=int(snk.params.ancestral_Ne),
        ratemap=ratemap, seed=int(snk.params.recapitate_seed)
    )
    
    ts = simulate_neutral_mutation(
        ts=ts, neutral_rate=float(snk.params.neutral_mu), seed=int(snk.params.neutral_seed)
    )
    
    ts = maf_threshold(ts, maf=float(snk.params.maf))

    ts.dump(snk.output[0])

if __name__ == '__main__':
    main()