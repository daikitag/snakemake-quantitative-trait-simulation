import tskit
import numpy as np
import pandas as pd

from snakemake.script import snakemake as snk


def update_result_dict(individual_ids, population, result_dict):
    for ind_id in individual_ids:
        result_dict["individual_id"].append(ind_id)
        result_dict["population"].append(population)
        result_dict["plink_id"].append(f"{population}_{ind_id}")


def create_df(yri_individuals, ceu_individuals):
    result_dict = {
        "individual_id": [],
        "plink_id": [],
        "population": [],
    }
    update_result_dict(individual_ids=yri_individuals,
                       population="YRI",
                       result_dict=result_dict)
    update_result_dict(individual_ids=ceu_individuals,
                       population="CEU",
                       result_dict=result_dict)
    
    return pd.DataFrame(result_dict)

def obtain_population_id(ts):
    for population in ts.populations():
        if population.metadata["name"] == "EXPLOSIVE_YRI":
            explosive_yri_id = population.id
        elif population.metadata["name"] == "EXPLOSIVE_CEU":
            explosive_ceu_id = population.id
            
    return explosive_yri_id, explosive_ceu_id

def obtain_individual_df(ts, yri_number, ceu_number, rng):
    explosive_yri_id, explosive_ceu_id = obtain_population_id(ts)
    
    yri_individuals = np.unique(ts.nodes_individual[ts.samples(explosive_yri_id)])
    ceu_individuals = np.unique(ts.nodes_individual[ts.samples(explosive_ceu_id)])
    
    yri_selected_individuals = np.sort(rng.choice(yri_individuals, yri_number, replace=False))
    ceu_selected_individuals = np.sort(rng.choice(ceu_individuals, ceu_number, replace=False))
    
    individual_id_df = create_df(yri_selected_individuals, ceu_selected_individuals)
    
    return individual_id_df

def main():    
    ts = tskit.load(snk.input[0])
    rng = np.random.default_rng(seed=int(snk.params.seed))

    individual_id_df = obtain_individual_df(
        ts = ts,
        yri_number = int(snk.params.yri_number),
        ceu_number = int(snk.params.ceu_number),
        rng = rng
    )
    
    individual_id_df.to_csv(snk.output[0], index=False)
    
if __name__ == '__main__':
    main()