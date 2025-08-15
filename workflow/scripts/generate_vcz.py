import tszip
import numpy as np
import pandas as pd

import tskit
import collections

import bio2zarr.tskit as ts2z

from snakemake.script import snakemake as snk


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
    ts = tszip.load(snk.input.ts)
    ts = maf_threshold(ts, maf=float(snk.params.maf))

    individual_id_df = pd.read_csv(snk.input.individual_id)

    model_mapping = ts.map_to_vcf_model(
        individuals=individual_id_df["individual_id"],
        individual_names=individual_id_df["plink_id"]
    )

    ts2z.convert(
        ts,
        vcz_path=snk.output[0],
        model_mapping=model_mapping,
        worker_processes=snk.threads
    )

    
if __name__ == '__main__':
    main()