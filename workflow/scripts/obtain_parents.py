import tskit
import pandas as pd
import numpy as np

from snakemake.script import snakemake as snk

def update_result_dict(ts, individual_ids, population, result_dict):
    for ind_id in individual_ids:
        individual = ts.individual(ind_id)
        result_dict["individual_id"].append(ind_id)
        result_dict["population"].append(population)
        result_dict["plink_id"].append(f"{population}_{ind_id}")
        result_dict["pedigree_p1"].append(individual.metadata["pedigree_p1"])
        result_dict["pedigree_p2"].append(individual.metadata["pedigree_p2"])

def obtain_population_id(ts):
    for population in ts.populations():
        if population.metadata["name"] == "EXPLOSIVE_YRI":
            explosive_yri_id = population.id
        elif population.metadata["name"] == "EXPLOSIVE_CEU":
            explosive_ceu_id = population.id
            
    return explosive_yri_id, explosive_ceu_id

def create_df(ts, yri_individuals, ceu_individuals):
    result_dict = {
        "individual_id": [],
        "plink_id": [],
        "population": [],
        "pedigree_p1": [],
        "pedigree_p2": [],
    }
    update_result_dict(ts, individual_ids=yri_individuals,
                       population="YRI",
                       result_dict=result_dict)
    update_result_dict(ts, individual_ids=ceu_individuals,
                       population="CEU",
                       result_dict=result_dict)
    
    return pd.DataFrame(result_dict)

def obtain_parent_df(ts):
    explosive_yri_id, explosive_ceu_id = obtain_population_id(ts)

    yri_individuals = np.unique(ts.nodes_individual[ts.samples(explosive_yri_id)])
    ceu_individuals = np.unique(ts.nodes_individual[ts.samples(explosive_ceu_id)])

    parent_df = create_df(ts, yri_individuals, ceu_individuals)

    return parent_df

def main():
    ts = tskit.load(snk.input[0])

    parent_df = obtain_parent_df(ts)

    parent_df.to_csv(snk.output[0], index=False)
    
if __name__ == '__main__':
    main()