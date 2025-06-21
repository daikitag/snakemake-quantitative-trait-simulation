import pandas as pd

from snakemake.script import snakemake as snk

def obtain_plink_phenotype(ts, phenotype_df):
    plink_phenotype_df = pd.DataFrame({
        "FID": [0 for _ in range(ts.num_individuals)],
        # PLINK doesn't like _
        "IID": [f"tsk_{i}indv" for i in range(ts.num_individuals)],
        "phenotype": phenotype_df.phenotype
    })
    
    return plink_phenotype_df

def main():    
    ts = tszip.load(snk.input.ts)
    
    phenotype_df = pd.read_csv(snk.input.phenotype_df)
    
    plink_phenotype = obtain_plink_phenotype(ts=ts, phenotype_df=phenotype_df)
    
    plink_phenotype.to_csv(snk.output[0], sep='\t', index=False)
    
if __name__ == '__main__':
    main()
