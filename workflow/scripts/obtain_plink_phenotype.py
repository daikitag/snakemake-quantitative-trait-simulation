import pandas as pd

from snakemake.script import snakemake as snk

def obtain_plink_phenotype(phenotype_df, individual_id_df):
    plink_phenotype_df = pd.DataFrame({
        "FID": individual_id_df.plink_id.tolist(),
        "IID": individual_id_df.plink_id.tolist(),
        "phenotype": pd.merge(individual_id_df, phenotype_df, on="individual_id").phenotype.to_list()
    })
    
    return plink_phenotype_df

def main():   
    phenotype_df = pd.read_csv(snk.input.phenotype_df)

    individual_id_df = pd.read_csv(snk.input.individual_id)
    
    plink_phenotype = obtain_plink_phenotype(
        phenotype_df=phenotype_df, individual_id_df=individual_id_df
    )
    
    plink_phenotype.to_csv(snk.output[0], sep='\t', index=False)
    
if __name__ == '__main__':
    main()
