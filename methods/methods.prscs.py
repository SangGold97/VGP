import os
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

from optparse import OptionParser
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore') 

def prs_cs(src='PRScs/PRScs.py',
            ref_dir=None,
            bim_prefix=None,
            sst_file=None,
            n_gwas=79550,
            chrom='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22',
            out_dir=None,
            phi=1e-2,
            seed=68):

    os.system(f'python {src} --ref_dir={ref_dir} \
                            --bim_prefix={bim_prefix} \
                            --sst_file={sst_file} \
                            --n_gwas={n_gwas} \
                            --chrom={chrom} \
                            --phi={phi} \
                            --seed={seed} \
                            --out_dir={out_dir}')

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-r", "--ref-dir", type=str, default='../reference_data/1kgp/ldblk_1kg_eas')
    parser.add_option("-b", "--bim-prefix", type=str, default=None)
    parser.add_option("-s", "--sst-file", type=str, default=None)
    parser.add_option("-n", "--n-gwas", type=str, default=None)
    parser.add_option("-i", "--phi", type=str, default=None)
    parser.add_option("-c", "--chrom", type=str, default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22')
    parser.add_option("-d", "--out-dir", type=str, default=None)
    
    (opt, args) = parser.parse_args()
    prs_cs(src='PRScs/PRScs.py',
            ref_dir=opt.ref_dir,
            bim_prefix=opt.bim_prefix,
            sst_file=opt.sst_file,
            n_gwas=opt.n_gwas,
            phi=opt.phi,
            chrom=opt.chrom,
            seed=68,
            out_dir=opt.out_dir)