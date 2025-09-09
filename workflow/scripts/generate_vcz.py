import numpy as np
import pandas as pd

import tskit
import collections

import bio2zarr.tskit as ts2z

#from snakemake.script import snakemake

def main():    
    ts = tskit.load(snakemake.input.ts)
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