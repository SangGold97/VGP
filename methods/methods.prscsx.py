import os
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

from optparse import OptionParser
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# PRS-CSx function
def prs_csx(src='PRScsx/PRScsx.py',
            ref_dir=None,
            bim_prefix=None,
            sst_file=None,
            n_gwas=None,
            pop=None,
            phi=1e-2,
            chrom=None,
            meta=True,
            seed=68,
            out_dir=None,
            out_name=None):

    os.system(f'python {src} --ref_dir={ref_dir} \
                            --bim_prefix={bim_prefix} \
                            --sst_file={sst_file} \
                            --n_gwas={n_gwas} \
                            --pop={pop} \
                            --phi={phi} \
                            --chrom={chrom} \
                            --meta={meta} \
                            --seed={seed} \
                            --out_dir={out_dir} \
                            --out_name={out_name}')
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-r", "--ref-dir", type=str, default=None)
    parser.add_option("-b", "--bim-prefix", type=str, default=None)
    parser.add_option("-s", "--sst-file", type=str, default=None)
    parser.add_option("-n", "--n-gwas", type=str, default=None)
    parser.add_option("-p", "--pop", type=str, default=None)
    parser.add_option("-i", "--phi", type=str, default=1e-2)
    parser.add_option("-c", "--chrom", type=str, default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22')
    parser.add_option("-d", "--out-dir", type=str, default=None)
    parser.add_option("-o", "--out-name", type=str, default=None)
    
    (opt, args) = parser.parse_args()
    prs_csx(src='PRScsx/PRScsx.py',
            ref_dir=opt.ref_dir,
            bim_prefix=opt.bim_prefix,
            sst_file=opt.sst_file,
            n_gwas=opt.n_gwas,
            pop=opt.pop,
            phi=opt.phi,
            meta=True,
            chrom=opt.chrom,
            seed=68,
            out_dir=opt.out_dir,
            out_name=opt.out_name)