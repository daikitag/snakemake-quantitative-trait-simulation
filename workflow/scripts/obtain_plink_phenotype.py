import pandas as pd

from snakemake.script import snakemake as snk

def obtain_plink_phenotype(phenotype_df, individual_id_df):
    plink_phenotype_df = pd.DataFrame({
        "FID": individual_id_df.plink_id.tolist(),
        "IID": individual_id_df.plink_id.tolist(),
        "phenotype": pd.merge(individual_id_df, phenotype_df, on="individual_id").phenotype.to_list()
    })
    
    return plink_phenotype_df

def merge_phenotype_df(phenotype_df_list):
    phenotype_df = pd.read_csv(phenotype_df_list[0])
    phenotype_df = phenotype_df.sort_values(by=["individual_id"])

    for i in range(1, len(phenotype_df_list)):
        add_df = pd.read_csv(phenotype_df_list[i])
        add_df = add_df.sort_values(by=["individual_id"])
        phenotype_df["phenotype"] += add_df["phenotype"]
    
    return phenotype_df

def main():
    phenotype_df = merge_phenotype_df(snk.input.phenotype_df_list)

    individual_id_df = pd.read_csv(snk.input.individual_id)
    
    plink_phenotype = obtain_plink_phenotype(
        phenotype_df=phenotype_df, individual_id_df=individual_id_df
    )
    
    plink_phenotype.to_csv(snk.output[0], sep='\t', index=False)
    
if __name__ == '__main__':
    main()
