# Script to calculate multispecies conservation L (D in the text). Requires 4nes.py to be run before.

import pandas as pd
import numpy as np
import math
from Bio import SeqIO
import csv
from Bio.Seq import Seq
from ete3 import Tree
from Bio import AlignIO
import matplotlib
from matplotlib import mlab
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.stats import entropy
import os
import sys

# to get alignment
def get_nt(tf, gene):

    code = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    syn_codons = ['A', 'G', 'P', 'T', 'V']

    all_ = pd.read_excel(tf + '/clustered.xlsx', index_col = [0])

    ids0 = list(all_.loc[gene]['bac'])

    outfile = tf + '/' + str(int(gene)) + '_prot_ali.fasta' # already there after 4nes.py
    align = AlignIO.read(outfile, "fasta")

    ids = [rec.id for rec in align]

    aa = pd.DataFrame(np.array([list(rec.seq) for rec in align], np.character))
    aa = aa.applymap(lambda x: str(x)[-2])
    aa = aa.set_index([pd.Index(ids)])
    #print ('got aa')

    dic_seq = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(tf + '/' + str(int(gene)) + '.fasta', 'fasta')}
    #print ('got dic_seq')

    nt_list = []
    for bac in ids:
        nt = []
        nucseq = dic_seq[bac]
        c = 0
        for l in aa.loc[bac]:
            if l != '-':
                if l != code[nucseq[:3]]:
                    print (gene, bac, 'MISMATCH!', l, code[nucseq[:3]])
                nt.append(nucseq[:3])
                nucseq = nucseq[3:]
            else:
                nt.append('-')
            c += 1
        nt_list.append(nt)
    nt = pd.DataFrame(np.array(nt_list))
    nt = nt.set_index([pd.Index(ids)])
    #print ('got nt')

    for i in range(aa.shape[1]):
        if len(set(aa[i])) > 1 or list(set(list(aa[i])) & set(syn_codons)) == []:
            del aa[i]
            del nt[i]
    nt = nt.applymap(lambda x: x[2])
    nt.columns = [i for i in range(aa.shape[1])]

    return [nt, ids]

# to get sites
def get_site_df(tf, gene, ids):
    dic_site0 = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(tf + '/' + str(int(gene)) + '_sites.fasta', 'fasta')}
    dic_site = {i: dic_site0[i] for i in ids}
    site_lst = [[bac] + list(dic_site[bac]) for bac in ids]
    site_df = pd.DataFrame(site_lst).set_index(0)
    site_df = site_df.applymap(lambda x: x.upper())
    return site_df

# to calculate multispecies conservation L (D in the text)
def get_L(tf, ids, site_df, nt, tree_file):

    def get_dist(node1, node2, tree):
        n1 = tree&node1
        n2 = tree&node2
        d = n1.get_distance(n2)
        return d

    t = Tree(tree_file)

    L_sites = []
    L_genes = []

    for n in range(1, len(dic_cons[tf])):
        for bac in ids:
            if dic_cons[tf][n] != '' and site_df.at[bac, n] != dic_cons[tf][n]:
                ref = bac
                noncons = site_df.at[bac, n]
                d = {key: value for key, value in zip([get_dist(ref, species, t) for species in ids], list(site_df[n]))}
                d_sorted = {key: d[key] for key in sorted(d.keys())}
                L = 0
                for k in d_sorted.keys():
                    if d_sorted[k] == d_sorted[0.0]:
                        L = k
                    if d_sorted[k] != d_sorted[0.0]:
                        break
                L_sites.append(L)
                for p in range(nt.shape[1]):
                    if nt.at[ref, p] == noncons:
                        d = {key: value for key, value in zip([get_dist(ref, species, t) for species in ids], list(nt[p]))}
                        d_sorted = {key: d[key] for key in sorted(d.keys())}
                        l = 0
                        for k in d_sorted.keys():
                            if d_sorted[k] == d_sorted[0.0]:
                                l = k
                            if d_sorted[k] != d_sorted[0.0]:
                                break
                        L_genes.append(l)

    return [L_sites, L_genes]

# plotting smooth distribution of L
def get_plot1(tf, L_sites, L_genes, min_spec):

    L_sites_dict = {key: value for key, value in zip([i for i in range(len(L_sites))], L_sites)}
    L_genes_dict = {key: value for key, value in zip([i for i in range(len(L_genes))], L_genes)}

    s = pd.Series(L_sites_dict)
    g = pd.Series(L_genes_dict)
    fig, ax = plt.subplots(figsize=(8, 4))

    ax = s.plot.kde(label='Observed', color='#fa8072')
    ax = g.plot.kde(label='Expected', color='skyblue')
    ax.grid(True)
    ax.legend(loc='right')
    ax.set_title('L distribution for ' + str(tf), fontsize=20)
    ax.legend(loc='upper right', fontsize=15)
    ax.set_xlabel('value of L', fontsize=15)
    ax.set_ylabel('Probability density', fontsize=15)
    ax.set_title('Min species '+str(min_spec+1), fontsize=16)
    if tf == '':
        ax.figure.savefig('pics/dist_smooth' + str(min_spec) + '.png', dpi=500, transparent=True)
    else:
        ax.figure.savefig(tf + '/dist_smooth' + str(min_spec) + '.png', dpi=500, transparent=True)

# plotting stepped distribution of L
def get_plot2(tf, L_sites, L_genes, min_spec, upper_range, bins):

    L_sites_dict = {key: value for key, value in zip([i for i in range(len(L_sites))], L_sites)}
    L_genes_dict = {key: value for key, value in zip([i for i in range(len(L_genes))], L_genes)}

    s = pd.Series(L_sites_dict)
    g = pd.Series(L_genes_dict)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax = g.plot.hist(bins = bins, histtype='stepfilled', alpha = 0.7, density=True, range=(0,upper_range), label='Expected', color='skyblue')
    ax = s.plot.hist(bins = bins, histtype='stepfilled', alpha = 0.4, density=True, range=(0,upper_range), label='Observed', color='#fa8072')
    ax.grid(True)
    ax.set_title('L distribution for ' + str(tf), fontsize=20)
    ax.legend(loc='upper right', fontsize=15)
    ax.set_xlabel('value of L', fontsize=15)
    ax.set_ylabel('Probability density', fontsize=15)
    ax.set_title('Min species '+str(min_spec+1), fontsize=16)
    if tf == '':
        ax.figure.savefig('pics/dist_step' + str(min_spec) + '.png', dpi=500, transparent=True)
    else:
        ax.figure.savefig(tf + '/dist_step' + str(min_spec) + '.png', dpi=500, transparent=True)

def KL_divergence(a, b):
    b = np.where(b == 0.0, 1e-6, b)
    return entropy(a, b)

# Big cycle for all TFs

tree_file = 'baci58_bacnames_rooted.nwk'

dic_cons = {'CcpA': {1: '', 2: 'T', 3: 'G', 4: '', 5: 'A', 6: 'A', 7: '', 8: 'C', 9: 'G', 10: '', 11: 'T', 12: 'T', 13: '', 14: 'C', 15: 'A', 16: ''},
            'CodY': {1: 'A', 2: 'A', 3: 'T', 4: 'T', 5: 'T', 6: 'T', 7: '', 8: '', 9: '', 10: 'A', 11: 'A', 12: 'A', 13: 'A', 14: 'T', 15: 'T'}, 
            'Fur': {1: 'A', 2: '', 3: 'T', 4: 'G', 5: 'A', 6: '', 7: 'A', 8: 'A', 9: 'T', 10: '', 11: 'A', 12: 'T', 13: 'T', 14: '', 15: 'T', 16: 'C', 17: 'A', 18: '', 19: 'T'},
            'LexA': {1: '', 2: '', 3: 'G', 4: 'A', 5: 'A', 6: 'C', 7: 'A', 8: 'T', 9: 'A', 10: 'T', 11: 'G', 12: 'T', 13: 'T', 14: 'C', 15: '', 16: ''},
            'PerR': {1: 'T', 2: 'T', 3: 'A', 4: 'T', 5: 'A', 6: 'A', 7: 'T', 8: '', 9: 'A', 10: 'T', 11: 'T', 12: 'A', 13: 'T', 14: 'A', 15: 'A'},
            'PurR': {1: 'A', 2: 'A', 3: 'A', 4: '', 5: '', 6: '', 7: 'G', 8: 'A', 9: 'A', 10: '', 11: 'A', 12: '', 13: 'T', 14: '', 15: '', 16: '', 17: '', 18: '', 19: '', 20: '', 21: '', 22: '', 23: '', 24: '', 25: '', 26: '', 27: '', 28: '', 29: '', 30: '', 31: '', 32: 'A', 33: '', 34: 'T', 35: '', 36: 'T', 37: 'T', 38: 'C', 39: '', 40: '', 41: '', 42: 'T', 43: 'T', 44: 'T'},
            'TnrA': {1: 'T', 2: 'G', 3: 'T', 4: '', 5: 'A', 6: '', 7: 'A', 8: '', 9: '', 10: '', 11: 'T', 12: '', 13: 'T', 14: '', 15: 'A', 16: 'C', 17: 'A'}}

L_all_sites = []
L_all_genes = []

Us = []
ps = []
KL = []

tfs = ['Fur', 'PerR', 'PurR', 'LexA', 'TnrA', 'CodY', 'CcpA']

min_spec = 9 # minimal alignment size

for tf in tfs:
    print ('TF', tf)

    f = pd.read_excel(tf + '/clustered.xlsx', index_col = [0])
    genes = list(set([i for i in list(f.index.dropna()) if list(f.index.dropna()).count(i) > min_spec]))

    print ('genes of this tf', genes)
    L_single_sites = []
    L_single_genes = []
    for gene in genes:
        print ('gene', gene)
        got_nt = get_nt(tf, gene)
        print ('got nt')
        nt = got_nt[0]
        ids = got_nt[1]
        site_df = get_site_df(tf, gene, ids)
        got_L = get_L(tf, ids, site_df, nt, tree_file)
        print ('got L')
        print ('number of dots for sites', len(got_L[0]))
        print ('number of dots for genes', len(got_L[1]))
        L_single_sites += got_L[0]
        L_single_genes += got_L[1]

    L_all_sites += L_single_sites
    L_all_genes += L_single_genes

    # for U-test
    U, p = mannwhitneyu(L_single_sites, L_single_genes)
    Us.append(U)
    ps.append(p)
    # for KL-distance
    #upper_range = np.max(L_single_sites + L_single_genes)
    upper_range = 2.8
    bins = 30
    ls = np.histogram(L_single_sites, bins=bins, range=(0,upper_range), density=True)[0]
    lg = np.histogram(L_single_genes, bins=bins, range=(0,upper_range), density=True)[0]
    kl = KL_divergence(ls, lg)
    KL.append(kl)

    # plot
    sns.set(font_scale=1.1)
    plt.style.use('seaborn-v0_8-whitegrid')
    try:
        get_plot1(tf, L_single_sites, L_single_genes, min_spec)
        get_plot2(tf, L_single_sites, L_single_genes, min_spec, upper_range, bins)
        # save L_single_sites, L_single_genes as lists
        L_single_sites_file = open(tf+'/L_single_sites_file.txt', 'w')
        L_single_sites_file.write(str(L_single_sites))
        L_single_sites_file.close()
        L_single_genes_file = open(tf+'/L_single_genes_file.txt', 'w')
        L_single_genes_file.write(str(L_single_genes))
        L_single_genes_file.close()
        print ('got plots')

    except TypeError:
        continue

final = {'tfs': tfs, 'U-test': Us, 'p-value': ps, 'KL-distance': KL}
final_df = pd.DataFrame(data=final)
final_df = final_df.set_index('tfs')
final_df.to_excel('final_df_L.xlsx')

get_plot1('', L_all_sites, L_all_genes, min_spec)
get_plot2('', L_all_sites, L_all_genes, min_spec, upper_range, bins)
# save L_all_sites, L_all_genes as lists
L_all_sites_file = open('L_all_sites_file.txt', 'w')
L_all_sites_file.write(str(L_all_sites))
L_all_sites_file.close()
L_all_genes_file = open('L_all_genes_file.txt', 'w')
L_all_genes_file.write(str(L_all_genes))
L_all_genes_file.close()


