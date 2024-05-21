# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 16:37:42 2023

@author: h
"""
#!/usr/bin/env python


import shutil
import argparse
from pyplink import PyPlink
import random
import subprocess
import numpy as np

def get_resampled_loci(number_markers):
    '''
    returns a list of loci positions resampled with replacement
    '''
    return sorted([random.randint(0, number_markers-1) for _ in range(number_markers)])

### parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('bfile', help='plink binary (bed, bim, fam) dataset prefix')
parser.add_argument('nboot', help='number of bootstrap replicates', type=int)
args = parser.parse_args()

bfile = r'F:\svdata\sheep_532_0515_geno025'
nboot = 1000
### perform plink bootstrap replicates
with PyPlink(bfile) as bed:
    bim = bed.get_bim() # returns pandas.DataFrame of bim file

    nmarkers = bed.get_nb_markers()
    nsamples = bed.get_nb_samples()
    print(f"### Loaded {nmarkers} markers and {nsamples} samples...")

    for rep in range(nboot):
        print(f"### Performing bootstrap replicate {rep+1} of {nboot}...")
        rep_list = get_resampled_loci(nmarkers) # gets a list of resampled markers

        rep_basename = f"rep{rep}"
        with PyPlink(rep_basename, "w") as outbed, open(f"{rep_basename}.bim", "w") as outbim:

            for marker_position in rep_list:                
                bed.seek(marker_position)
                marker, genotypes = next(bed)
                # write marker info to outbim
                marker_info = bim.loc[marker]
                marker_line = '\t'.join(map(str, marker_info.tolist())) + '\n'
                outbim.write(marker_line)
                # write genotypes to outbed
                outbed.write_genotypes(genotypes)

        # create copy of fam file with rep_basename.fam
        # 使用 Python 的 shutil.copyfile 来复制文件
        fam_original = f"{bfile}.fam"
        fam_copy = f"{rep_basename}.fam"
        shutil.copyfile(fam_original, fam_copy)
        # call plink to compute distance matrix on rep dataset
        plink_cmd_line = f"D:\\下载\\plink_win64_20231211\\plink --bfile {rep_basename} --distance square 1-ibs flat-missing --out {rep_basename} --chr-set 27 --allow-extra-chr"
        subprocess.run(plink_cmd_line, shell=True, stdout=subprocess.DEVNULL)
        # clean created files, leaving only matrices (.mdist)
        #clean_cmd_line = f"rm {rep_basename}.bed {rep_basename}.bim {rep_basename}.fam {rep_basename}.nosex {rep_basename}.log {rep_basename}.mdist.id"
        #subprocess.run(clean_cmd_line.split())

### concatenate and format all generated matrices for input to PHYLIP

# get sample names from .fam file
samples = []
with open(f"{bfile}.fam", "r") as fam:
    for line in fam:
        sample_id = line.strip().split()[1]
        samples.append(sample_id)

with open("infile_sheep1", "w") as outf:

    for rep in range(nboot):
        if rep == 989 or rep == 989 or rep == 330 or rep == 329 or rep == 328 or rep == 327 or rep == 326 or rep == 325 or rep == 324 or rep == 322 or rep == 249 or rep == 248 or rep == 247 or rep == 242:
            continue
        # load each replicate matrix as a numpy array
        rep_file = f"rep{rep}.mdist"
        rep_matrix = np.loadtxt(rep_file)

        # write matrix in PHYLIP format
        outf.write(f"    {len(samples)}\n")
        for i, sample in enumerate(samples):
            matrix_line = '  '.join(f"{x:0.6f}" for x in rep_matrix[i])
            out_line = f"{sample.ljust(12)}{matrix_line}\n"
            outf.write(out_line)
        outf.write("\n")


infile = pd.read_csv(r'D:\scpy\infile', sep='\t',header=None)
infile_new = pd.read_csv(r'C:\Users\h\infile_0516_capra', sep='\t',header=None)

