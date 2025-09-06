import numpy as np
import pandas as pd

import tskit
import collections

import bio2zarr.tskit as ts2z

#from snakemake.script import snakemake


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

def subset_tree_seq(ts, selected_individuals):
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate((selected_nodes, ts.individual(individual).nodes))

    return ts.simplify(selected_nodes, filter_individuals=False)

def main():    
    ts = tskit.load(snakemake.input.ts)
    individual_id_df = pd.read_csv(snakemake.input.individual_id)
    ts = subset_tree_seq(ts, individual_id_df["individual_id"])    
    ts = maf_threshold(ts, maf=float(snakemake.params.maf))
   
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