import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit

from snakemake.script import snakemake as snk

def convert_allele(ts, seed):
    ts = pyslim.generate_nucleotides(ts, seed=seed)
    ts = pyslim.convert_alleles(ts)
    
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


def subset_tree_seq(ts, selected_individuals):
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate((selected_nodes, ts.individual(individual).nodes))

    return ts.simplify(selected_nodes, filter_individuals=False)

def main():
    
    ts = tskit.load(snk.input.ts)
    
    # Nucleotide instead of SLiM allele
    ts = convert_allele(ts, seed=int(snk.params.chromosome))
    ratemap = obtain_ratemap(snk.params.genetic_map)
    
    # Recapitate
    ts = recapitate(
        ts=ts, ancestral_Ne=int(snk.params.ancestral_Ne),
        ratemap=ratemap,
        seed=int(snk.params.recapitate_seed)+int(snk.params.chromosome)
    )

    individual_id_df = pd.read_csv(snk.input.individual_id)
    ts = subset_tree_seq(ts, individual_id_df["individual_id"])
    

    ts.dump(snk.output[0])

if __name__ == '__main__':
    main()