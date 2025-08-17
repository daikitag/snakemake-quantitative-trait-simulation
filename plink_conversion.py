import sys
from vcztools.plink import write_plink

write_plink(vcz_path=sys.argv[1], out=sys.argv[2])