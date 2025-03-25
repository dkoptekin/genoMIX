#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

"""
@author: Dilek Koptekin
"""

import yaml
import numpy as np
import pandas as pd
import itertools
import math
import os
import sys

# TODO add merge option or keep_pops
# TODO add option to just using overlap SNPs

def readInd(indPath):
    return np.loadtxt(indPath,dtype=str).transpose()

def readSnp(snpPath):
    return pd.read_csv(snpPath, delim_whitespace=True, header=None, names = ["snpID", "chr", "cM", "pos", "allele1", "allele2"])

def readUnpackedgeno(genoPath):
    return pd.read_table(genoPath, header=None, names=["genotypes"], dtype=str, iterator=True, low_memory=False)

def readUnpackedgenoChunks(genoChunk, nInd):
    g = np.array(genoChunk.genotypes, dtype=f'<U{nInd}')
    return g.view('int32').reshape((g.size, -1)) - 48

def readPackedgeno(genoPath):
    genoTmp = open(genoPath, "rb")
    header = genoTmp.read(20)
    nInd, nSnp = [int(x) for x in header.split()[1:3]]
    rLen = max(48, int(np.ceil(nInd * 2 / 8)))
    genoTmp.close()
    genoFile = np.fromfile(genoPath, dtype='uint8')[rLen:]  # without header
    genoFile.shape = (nSnp, rLen)
    return genoFile, nInd, nSnp

def readPackedgenoChunks(genoChunk, nInd):
    g = np.unpackbits(genoChunk, axis=1)[:, :(2 * nInd)]
    g = 2 * g[:, ::2] + g[:, 1::2]
    g[g == 3] = 9
    return g

def loadYaml(yamlFile):
    params = yaml.load(open(yamlFile), Loader=yaml.FullLoader)
    return params

def defineChunks(snpList, chunkSize = 1000):
    nChunks = len(snpList) // chunkSize + 1
    return np.array_split(snpList, nChunks)

def decodeParams(parameterDict, nWindow):
    source_pops, sources, rates, samplings  = [], [], [], []
    if '2-way' in parameterDict.keys():
        pop2s = parameterDict.get('2-way')[0].get('sources')
        rates2s = parameterDict.get('2-way')[1].get('rates')
        sampling2s = [np.random.choice(a=pop2s, size=nWindow, p=rates2s[r2]) for r2 in range(len(rates2s))]
        source2s = [pop2s] * len(rates2s)
        source_pops.append(pop2s), rates.append(rates2s), samplings.append(sampling2s), sources.append(source2s)
    if '3-way' in parameterDict.keys():
        pop3s = parameterDict.get('3-way')[0].get('sources')
        rates3s = parameterDict.get('3-way')[1].get('rates')
        sampling3s = [np.random.choice(a=pop3s, size=nWindow, p=rates3s[r3]) for r3 in range(len(rates3s))]
        source3s = [pop3s] * len(rates3s)
        source_pops.append(pop3s), rates.append(rates3s), samplings.append(sampling3s), sources.append(source3s)
    if '4-way' in parameterDict.keys():
        pop4s = parameterDict.get('4-way')[0].get('sources')
        rates4s = parameterDict.get('4-way')[1].get('rates')
        sampling4s = [np.random.choice(a=pop4s, size=nWindow, p=rates4s[r4]) for r4 in range(len(rates4s))]
        source4s = [pop4s] * len(rates4s)
        source_pops.append(pop4s), rates.append(rates4s), samplings.append(sampling4s), sources.append(source4s)
    return list(set(itertools.chain(*source_pops))), list(itertools.chain(*sources)), list(itertools.chain(*rates)), list(itertools.chain(*samplings))

def mixturePops(genoChunk, sourcePop, indFile, popFile):
    scol = popFile[sourcePop][np.random.choice(range(len(popFile[sourcePop])))]
    icol = int(np.flatnonzero(indFile[0] == scol))
    return genoChunk[:,icol]

param_dict = loadYaml(sys.argv[1])

param_dict_keys = param_dict.keys()
data_prefix = param_dict['data_prefix']
data_type = param_dict['data_type']
chunk_size = param_dict['chunk_size']

geno_path=os.path.join(data_prefix + ".geno")
snp_path=os.path.join(data_prefix + ".snp")
ind_path=os.path.join(data_prefix + ".ind")

if os.path.exists(geno_path) and os.path.exists(snp_path) and os.path.exists(ind_path):
    print(f"The input Eigenstrat files ({data_type}): \n"
          f" {geno_path} \n"
          f" {snp_path} \n"
          f" {ind_path}")
    ind_file = readInd(ind_path)
    if data_type == 'unpacked':
        geno_file = readUnpackedgeno(geno_path)
        n_ind = ind_file.shape[1]
        n_snp = sum(1 for row in open(snp_path, 'r'))
        n_chunks = math.ceil(n_snp / chunk_size)
    elif data_type == 'packed':
        geno_file, n_ind, n_snp = readPackedgeno(geno_path)
        n_chunks = math.ceil(n_snp / chunk_size)
else:
    print("Invalid datatype! \n"
          "geno/snp/ind file should be in same directory and start with same prefix and \n"
          "data type should be either 'packed' or 'unpacked'")
    sys.exit()

windows = defineChunks(range(n_snp), chunk_size)
n_window = len(windows)

source_pops, admix_sources, admix_rates, sampling_list = decodeParams(param_dict, n_window)

admix_indv= [f'adm{len(source)}pops_{"_".join(str(int(r*100)) for r in rate)}' for source, rate, sampling in zip(admix_sources, admix_rates, sampling_list)]
admix_geno = pd.DataFrame(index=range(n_snp), columns=admix_indv)
indices = np.flatnonzero(np.isin(ind_file[2], source_pops))
subset_ind_file = ind_file[:,indices]
pop_file = {key:subset_ind_file[0][np.where(subset_ind_file[2]==key)].tolist() for key in subset_ind_file[2]}
print("genoMIX running \n")

for chunk in range(n_window):
    if data_type == "unpacked":
        geno_chunk = readUnpackedgenoChunks(geno_file.get_chunk(len(windows[chunk])))
    elif data_type == "packed":
        geno_chunk = readPackedgenoChunks(geno_file[windows[chunk]], n_ind)
    geno_chunk_subset = geno_chunk[:, indices]
    for sources, rates, samplings in zip(admix_sources, admix_rates, sampling_list):
        indv_id = f'adm{len(sources)}pops_{"_".join(str(int(r * 100)) for r in rates)}'
        s_pop = samplings[chunk]
        admix_chunk = mixturePops(geno_chunk_subset, s_pop, subset_ind_file, pop_file)
        admix_geno.loc[windows[chunk], indv_id] = admix_chunk
    print(end="\r%6.2f %%" % (chunk / (n_window - 1) * 100))

#admix_geno.shape
#admix_geno.isnull().values.any()

# write outputs
if ('output_prefix' in param_dict_keys) and (param_dict['output_prefix'] != None):
    output_prefix=param_dict['output_prefix']
else:
    output_prefix = "toyadmix"

out_geno_path=os.path.join(output_prefix + ".geno")
out_snp_path=os.path.join(output_prefix + ".snp")
out_ind_path=os.path.join(output_prefix + ".ind")

print(f"The output Eigenstrat files: \n"
      f" {out_geno_path} \n"
      f" {out_snp_path} \n"
      f" {out_ind_path} \n"
      f" will be written as unpacked")

# output geno file
with open(out_geno_path, 'w') as out_geno:
    for row in admix_geno.values:
        out_geno.write("".join(map(str, row)) +"\n")

#output ind file
sampleIDs = admix_geno.columns
nindv= len(sampleIDs)
admix_ind = pd.DataFrame({'col1': sampleIDs,
                          'col2': np.repeat("U", nindv),
                          'col3':sampleIDs})
admix_ind.to_csv(out_ind_path, header=False, sep='\t', index=False, encoding='utf-8')

# output snp file
out_snp = os.system(f"cp {snp_path} {out_snp_path}")

if out_snp == 0:
    print("Done")
else:
    print("Something went wrong")
