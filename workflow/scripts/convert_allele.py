import pyslim
import tszip

from snakemake.script import snakemake as snk

def convert_allele(ts):
    ts = pyslim.generate_nucleotides(ts)
    ts = pyslim.convert_alleles(ts)
    
    return ts

def main():    
    ts = tszip.load(snk.input[0])
    
    ts = convert_allele(ts)
    
    tszip.compress(ts, snk.output[0])
    
if __name__ == '__main__':
    main()
