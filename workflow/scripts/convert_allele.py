import pyslim
import tskit

from snakemake.script import snakemake as snk

def convert_allele(ts):
    ts = pyslim.generate_nucleotides(ts)
    ts = pyslim.convert_alleles(ts)
    
    return ts

def main():    
    ts = tskit.load(snk.input[ts])
    
    ts = convert_allele(ts)
    
    ts.dump(snk.output[0])
    
if __name__ == '__main__':
    main()
