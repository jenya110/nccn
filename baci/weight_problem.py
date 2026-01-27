# Investigating how selection acts on the energy (weight) of binding between TF and TFBS. Script with the shuffling procedure, and plotting of SD distributions and rank distributions.

import pandas as pd
import numpy as np
import math
import csv
import matplotlib
from matplotlib import mlab
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.stats import entropy


plt.style.use('seaborn-whitegrid')

dic_cons = {'CcpA': {1: '', 2: 'T', 3: 'G', 4: '', 5: 'A', 6: 'A', 7: '', 8: 'C', 9: 'G', 10: '', 11: 'T', 12: 'T', 13: '', 14: 'C', 15: 'A', 16: ''},
            'CodY': {1: 'A', 2: 'A', 3: 'T', 4: 'T', 5: 'T', 6: 'T', 7: '', 8: '', 9: '', 10: 'A', 11: 'A', 12: 'A', 13: 'A', 14: 'T', 15: 'T'}, 
            'Fur': {1: 'A', 2: '', 3: 'T', 4: 'G', 5: 'A', 6: '', 7: 'A', 8: 'A', 9: 'T', 10: '', 11: 'A', 12: 'T', 13: 'T', 14: '', 15: 'T', 16: 'C', 17: 'A', 18: '', 19: 'T'},
            'LexA': {1: '', 2: '', 3: 'G', 4: 'A', 5: 'A', 6: 'C', 7: 'A', 8: 'T', 9: 'A', 10: 'T', 11: 'G', 12: 'T', 13: 'T', 14: 'C', 15: '', 16: ''},
            'PerR': {1: 'T', 2: 'T', 3: 'A', 4: 'T', 5: 'A', 6: 'A', 7: 'T', 8: '', 9: 'A', 10: 'T', 11: 'T', 12: 'A', 13: 'T', 14: 'A', 15: 'A'},
            'PurR': {1: 'A', 2: 'A', 3: 'A', 4: '', 5: '', 6: '', 7: 'G', 8: 'A', 9: 'A', 10: '', 11: 'A', 12: '', 13: 'T', 14: '', 15: '', 16: '', 17: '', 18: '', 19: '', 20: '', 21: '', 22: '', 23: '', 24: '', 25: '', 26: '', 27: '', 28: '', 29: '', 30: '', 31: '', 32: 'A', 33: '', 34: 'T', 35: '', 36: 'T', 37: 'T', 38: 'C', 39: '', 40: '', 41: '', 42: 'T', 43: 'T', 44: 'T'},
            'TnrA': {1: 'T', 2: 'G', 3: 'T', 4: '', 5: 'A', 6: '', 7: 'A', 8: '', 9: '', 10: '', 11: 'T', 12: '', 13: 'T', 14: '', 15: 'A', 16: 'C', 17: 'A'}}

def get_E_dist(df, pwm):
    energies = []
    for i in df.index:
        site = ''.join(list(df.loc[i]))
        e = E(site, pwm)
        energies.append(e)
    return energies
    
def E(site, pwm):
    site_E = []
    for i in range(len(site)):
        site_E.append(pwm.at[site[i], i])
    #print (site_E)
    E = round(np.sum(site_E), 2)
    return E

def KL_divergence(a, b):
    b = np.where(b == 0.0, 1e-6, b)
    return entropy(a, b)

# Getting PWMs

tfs = ['Fur', 'PerR', 'PurR', 'LexA', 'TnrA', 'CodY', 'CcpA']

for tf in tfs:
    print(tf)

    f = pd.read_excel(tf + '/output.xlsx', index_col = [0])
    f = f.dropna(how='all')

    # wrighting ALL the sites into .fasta (including paralogs, so taking output.xlsx rather than clustered.xlsx)
    sites_file = open('weight_problem/' + tf + '/reg_sites.fasta', 'w')
    for i in list(f.site):
        sites_file.write('>\n' + i + '\n')
    sites_file.close()

    # making pwm
    length = len(i)

    sites_file = open('weight_problem/'+tf+'/reg_sites.fasta', 'r')
    l = [list(i.strip().lower()) for i in sites_file.readlines() if '>' not in i]
    sites = pd.DataFrame(l)

    pwm = pd.DataFrame(columns=[j for j in range(length)], index=['a', 'c', 'g', 't'])
    for i in range(sites.shape[1]):
        pwm[i] = [list(sites[i]).count('a')+1, list(sites[i]).count('c')+1, list(sites[i]).count('g')+1, list(sites[i]).count('t')+1]

    pwm1 = (pwm/(sites.shape[0]+4)).applymap(lambda x: math.log2(x))
    pwm1 = pwm1*-1
    for i in pwm1.columns:
        pwm1[i] = pwm1[i] - np.min(pwm1[i])

    pwm1.to_csv('weight_problem/'+tf+'/pwm_new')
    
    
# Getting SD distributions and rank distributions for all TFs

min_spec = 9 # minimal number of species in the alignments we consider

tfs = ['Fur', 'PerR', 'PurR', 'LexA', 'TnrA', 'CodY', 'CcpA']

tfs_sds = {}
tfs_mean_sds = {}
tfs_ranks = {}
tfs_energy_dists = {}
Us = {}
ps = {}
KL = {}

for tf in tfs:
    print('tf', tf)
    
    # get pwm
    pwm = pd.read_csv('weight_problem/'+tf+'/pwm_new')
    pwm = pwm.set_index('Unnamed: 0')
    pwm.index.name = None
    pwm.index = ['a', 'c', 'g', 't']
    pwm.columns = [0]+list(dic_cons[tf].keys())[:-1]
    
    # get sites
    f = pd.read_excel(tf + '/clustered.xlsx', index_col = [0])  
    genes = list(set([i for i in list(f.index.dropna()) if list(f.index.dropna()).count(i) > min_spec]))
    
    print(genes)
    
    f = f.reset_index()
    f = f[(f['po'].isin(genes))]

    f1 = f.drop(['gene_name', 'locus_tag', 'protseq', 'nucseq'], axis=1)
    f1['num'] = [i for i in range(f1.shape[0])]
    f1 = f1.set_index(['num'])
    f2 = pd.DataFrame([list(i) for i in list(f['site'])]) # сайты по буковкам
    f2['num'] = [i for i in range(f2.shape[0])]
    f2 = f2.set_index(['num'])
    
    df = pd.concat([f1,f2],axis=1,ignore_index=True)
    df.to_csv('weight_problem/' + tf + '/sites.csv' + str(min_spec))
    
    # consider genes separately
    
    mean_energies = []
    all_energies = []
    true_sds = []
    all_shuffled_sds = []
    ranks = []
        
    for gene in genes:
        df_gene = df[df[0] == gene].drop([0,1,2,3], axis=1)
        energies = get_E_dist(df_gene,pwm)
        all_energies.extend(energies) # adding E
        
        mean_energy = round(np.mean(energies), 2)
        mean_energies.append(mean_energy)
        true_sd = round(np.std(energies), 2)
        true_sds.append(true_sd)
        
        sds_shuffling = []
        
        for i in range(100): # shuffling 100 times
            df_gene_shuffled = df_gene.apply(np.random.permutation, axis=0)
            energies = get_E_dist(df_gene_shuffled,pwm)
            sd = round(np.std(energies), 2)
            sds_shuffling.append(sd)
        
        sds_shuffling_sorted = sorted(sds_shuffling)
        same_sds = []
        for e, j in enumerate(sds_shuffling_sorted):
            if j == true_sd:
                same_sds.append(e)
            if j > true_sd:
                break
        if same_sds == []:
            rank = e+1 # because it enumerates from 0
        else:
            rank = np.mean(same_sds)+1
        
        ranks.append(rank)
        all_shuffled_sds.extend(sds_shuffling_sorted) # добавить sds

    tfs_sds[tf] = round(np.std(mean_energies), 2)
    tfs_mean_sds[tf] = round(np.mean(true_sds), 2)
    
    tfs_ranks[tf] = ranks
    
    # Mann-Whitney test and p-value
    U, p = mannwhitneyu(true_sds, all_shuffled_sds, alternative='two-sided')
    Us[tf] = U
    ps[tf] = p
    
    # KL distance
    bins = 30
    upper_range = np.max(true_sds+all_shuffled_sds)
    ls = np.histogram(true_sds, bins=bins, range=(0,upper_range), density=True)[0]
    lg = np.histogram(all_shuffled_sds, bins=bins, range=(0,upper_range), density=True)[0]  
    kl = KL_divergence(ls, lg)
    KL[tf] = kl
    
    # plotting SDs
    plt.rcParams['font.size'] = '14'
    plt.figure(figsize=(4, 5))

    plt.hist(true_sds, bins=7, histtype='stepfilled', density=True, alpha=0.5, color='tomato', ec='k', label='real SD')
    plt.hist(all_shuffled_sds, bins=7, histtype='stepfilled', density=True, alpha=0.5, color='teal', ec='k', label='shuffled SD')
    plt.xlabel('SD')
    plt.ylabel('density')
    plt.legend()
    plt.title(f'{tf}, {len(genes)} genes')
    plt.tight_layout()
    plt.savefig(f'weight_problem/{tf}/sd1_sd2_9_new.png', dpi=400, transparent=True)
    plt.savefig(f'weight_problem/pics/{tf}_sd1_sd2_9_new.png', dpi=400, transparent=True)
    
    # ploting ranks
    plt.rcParams['font.size'] = '14'
    plt.figure(figsize=(2, 2.5))

    plt.hist(ranks, bins=7, histtype='stepfilled', alpha=0.5, color='violet', ec='k')
    plt.xlabel('ranks')
    plt.ylabel('count')
    plt.xlim(0,100)
    plt.xticks(np.arange(0, 101, step=50))
    plt.title(f'{tf}, {len(genes)} genes')
    plt.tight_layout()
    plt.savefig(f'weight_problem/{tf}/ranks_distribution_9_new.png', dpi=400, transparent=True)
    plt.savefig(f'weight_problem/pics/{tf}_ranks_distribution_9_new.png', dpi=400, transparent=True)
    

final = [Us, ps, KL]
final_df = pd.DataFrame.from_dict(final).T
final_df.columns = ['U-test', 'p-value', 'KL-distance']
final_df.to_excel('weight_problem/final_df_new.xlsx')

