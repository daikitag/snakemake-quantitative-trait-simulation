import numpy as np
import pandas as pd

import tskit
import pyslim

from snakemake.script import snakemake as snk

def convert_allele(ts):
    ts = pyslim.generate_nucleotides(ts)
    ts = pyslim.convert_alleles(ts)
    
    return ts

def subset_tree_seq(ts, selected_individuals):
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate((selected_nodes, ts.individual(individual).nodes))

    return ts.simplify(selected_nodes, filter_individuals=False)

def main():    
    ts = tskit.load(snk.input.ts)
    individual_id_df = pd.read_csv(snk.input.individual_id)
    ts = subset_tree_seq(ts, individual_id_df["individual_id"])
    ts = convert_allele(ts)
    ts.dump(snk.output[0])

    
if __name__ == '__main__':
    main()