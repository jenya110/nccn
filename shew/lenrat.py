# Script to calculate length ratio R. Requires 4nes.py to be run before.

import pandas as pd
import numpy as np
import math
from ete3 import Tree
import seaborn as sns
from Bio import SeqIO
import csv
from Bio.Seq import Seq
from Bio import AlignIO
import sys
import matplotlib
from matplotlib import mlab
import matplotlib.pyplot as plt
import seaborn as sns
import os

# correct bacterial names
dic_patric0 = {'Shewanella_sp._SNU_WT4_strain_SNU_WT1_chromosome': '2590015.3',
               'Shewanella_dokdonensis_strain_DSM_23626_chromosome': '712036.3',
               'Shewanella_sp._FJAT-51800_chromosome': '2814294.3',
               'Shewanella_aestuarii_strain_PN3F2_chromosome_PN3F2': '1028752.3',
               'Shewanella_maritima_strain_D4-2_chromosome': '2520507.3',
               'Shewanella_benthica_strain_DB21MT-2_genome_assembly': '43661.3',
               'Shewanella_donghaensis_strain_LT17_chromosome': '238836.3',
               'Shewanella_litorisediminis_strain_SMK1-12_chromosome': '1173586.3',
               'Shewanella_polaris_strain_SM1901_chromosome': '2588449.3',
               'Shewanella_khirikhana_strain_TH2012_chromosome': '1965282.3',
               'Shewanella_sp._MEBiC00475_chromosome': '2575361.3',
               'Shewanella_indica_strain_Colony474_chromosome': '768528.7',
               'Shewanella_livingstonensis_strain_LMG_19866_chromosome': '150120.4',
               'Shewanella_amazonensis_SB2B': '326297.10',
               'Shewanella_sp._KX20019_chromosome': '2803864.3',
               'Shewanella_denitrificans_OS217': '318161.16',
               'Shewanella_loihica_PV-4': '323850.11',
               'Shewanella_marisflavi_strain_EP1_chromosome': '260364.3',
               'Shewanella_sp._M2_chromosome': '2487742.3',
               'Shewanella_sp._WPAGA9_chromosome': '2775245.3',
               'Shewanella_sp._LZH-2_chromosome': '2806008.3',
               'Shewanella_sp._Arc9-LZ_chromosome': '2698686.3',
               'Shewanella_sp._Scap07_chromosome': '2589987.3',
               'Shewanella_decolorationis_strain_sesselensis_chromosome': '256839.7',
               'Shewanella_chilikensis_strain_DC57_chromosome': '558541.8',
               'Shewanella_sp._MR-4': '60480.19',
               'Shewanella_sp._MR-7': '60481.13',
               'Shewanella_violacea_DSS12': '637905.5',
               'Shewanella_frigidimarina_NCIMB_400': '318167.14',
               'Shewanella_japonica_strain_KCTC_22435_chromosome': '93973.3',
               'Shewanella_sp._W3-18-1': '351745.9',
               'Shewanella_pealeana_ATCC_700345': '398579.7',
               'Shewanella_baltica_strain_CW2_chromosome': '62322.6',
               'Shewanella_algae_strain_RQs-106_chromosome': '38313.77',
               'Shewanella_oneidensis_MR-1': '211586.12',
               'Shewanella_sp._Pdp11_chromosome': '2059264.4',
               'Shewanella_halifaxensis_HAW-EB4': '458817.8',
               'Shewanella_sp._ANA-3_chromosome_1': '94122.6',
               'Shewanella_bicestrii_strain_JAB-1_chromosome': '2018305.3',
               'Shewanella_sediminis_HAW-EB3': '425104.7',
               'Shewanella_sp._YLB-06_chromosome': '2593655.3',
               'Shewanella_piezotolerans_WP3': '225849.4',
               'Shewanella_sp._LC6_chromosome': '2589790.3',
               'Shewanella_putrefaciens_strain_SA70_chromosome': '24.6',
               'Shewanella_sp._WE21_chromosome': '2029986.5',
               'Shewanella_woodyi_ATCC_51908': '392500.6',
               'Shewanella_psychrophila_strain_WP2_chromosome': '225848.8'}

# to get protein sequence
def get_nt(tf, gene, all_):
    
    ids0 = list(all_.loc[gene]['bac']) # ids with spaces instead of _
    ids = ['_'.join(bac.split(' ')) for bac in ids0]
    
    gene_fasta_file = open(tf + '/' + str(int(gene)) + '.fasta', 'w') # nucseq
    site_fasta_file = open(tf + '/' + str(int(gene)) + '_sites.fasta', 'w') # site
    prot_fasta_file = open(tf + '/' + str(int(gene)) + '_prot.fasta', 'w') # protseq
    
    for bac in ids0:
        name = '_'.join(bac.split(' '))
        
        site = all_.loc[(all_['bac'] == bac)].at[gene, 'site']
        nucseq = all_.loc[(all_['bac'] == bac)].at[gene, 'nucseq']
        protseq = all_.loc[(all_['bac'] == bac)].at[gene, 'protseq']
        gene_fasta_file.write('>' + name + '\n') 
        gene_fasta_file.write(nucseq + '\n')
        prot_fasta_file.write('>' + name + '\n') 
        prot_fasta_file.write(protseq + '\n')
        site_fasta_file.write('>' + name + '\n')
        site_fasta_file.write(site + '\n')
        
    gene_fasta_file.close()
    site_fasta_file.close()
    prot_fasta_file.close()

    # aminoacid alignment
    infile = tf + '/' + str(int(gene)) + '_prot.fasta'
    outfile = tf + '/' + str(int(gene)) + '_prot_ali.fasta'
    bashCommand = f'muscle -align {infile} -output {outfile} -quiet'
    os.system(bashCommand)
    dict_prot = {}
    for protseq in SeqIO.parse(str(tf)+'/'+str(gene)+'_prot_ali.fasta', 'fasta'):
        dict_prot[protseq.id] = str(protseq.seq)
    outright = open(str(tf)+'/'+str(gene)+'_prot_ali.fasta', 'w')
    for i in ids:
        outright.write('>'+i+'\n'+dict_prot[i]+'\n')
    outright.close()


# to make .nuc file
def get_nuc(tf, gene):

    # count sequences
    file = open(str(tf)+'/'+str(gene)+'_prot_ali.fasta', 'r')
    lines = ''.join(file.readlines())
    len_seq = len(''.join(lines.split('>')[1].split('\n')[1:]))
    num_seq = lines.count('>')
    file.close()

    bacs = []
    for l in lines.split('>')[1:]:
        bacs.append(dic_patric0[l.split()[0]])
    long_bac = max([len(bac) for bac in bacs]) + 3
    
    codon_file = open(str(tf)+'/'+str(gene)+'.nuc', 'w')
    codon_file.write('  '+str(num_seq)+'   '+str(len_seq*3)+'\n')
    for protseq, nucseq in zip(SeqIO.parse(str(tf)+'/'+str(gene)+'_prot_ali.fasta', 'fasta'), SeqIO.parse(str(tf)+'/'+str(gene)+'.fasta', 'fasta')):
        codon_file.write(dic_patric0[protseq.id]+' '*(long_bac-len(dic_patric0[protseq.id])))
        codonseq = ''
        c = 0 # gaps counter
        for i in range(len(str(protseq.seq))):
            if protseq[i] != '-':
                codonseq += str(nucseq.seq)[(i*3-c*3):(i*3-c*3)+3]+' '
            if protseq[i] == '-':
                codonseq += '--- '
                c+=1
        codon_file.write(codonseq+'\n')
    codon_file.close()
    return bacs

# to rename bacteria from _sites.fasta to patric codes and transforms to nuc
def get_nuc_sites(tf, gene):
    
    sfile = open(tf+'/'+str(gene)+'_sites.fasta', 'r')
    a = ''.join(sfile.readlines())
    sfile.close()
        
    num_seq = a.count('>')
    len_seq = len(a.split('>')[1].split('\n')[-2])
    for k in dic_patric0.keys():
        a = a.replace(k, dic_patric0[k])

    bacs = []
    for l in a.split('>')[1:]:
        bacs.append(l.split()[0])
    long_bac = max([len(bac) for bac in bacs]) + 3

    codon_file = open(tf+'/'+str(gene)+'_sites.nuc', 'w')
    codon_file.write('  '+str(num_seq)+'   '+str(len_seq)+'\n')
    for i in a.split('>')[1:]:
        codon_file.write(i.split()[0]+' '*(long_bac-len(i.split()[0]))+i.split()[1]+'\n')
    codon_file.close()

# to rewrite ctl file for a new run of baseml/codeml
def rewrite_ctl(tf, gene, sites):
    if sites == False:
        inf = open('codeml.ctl', 'r')
        ctl = inf.readlines()
        inf.close()
        seqfile_line = ctl[0].split('= ')[0]+'= '+tf+'/'+str(gene)+'.'+ctl[0].split('.')[1]        
        treefile_line = ctl[1].split('= ')[0]+'= '+str(gene)+'tree'+ctl[1].split('tree')[-1]
        ctl_new = [seqfile_line, treefile_line]
        ctl_new += ctl[2:]
        
        outf = open('codeml.ctl', 'w')
        for line in ctl_new:
            outf.write(line)
        outf.close()
    if sites == True:
        inf = open('baseml.ctl', 'r')
        ctl = inf.readlines()
        inf.close()
        seqfile_line = ctl[0].split('= ')[0]+'= '+tf+'/'+str(gene)+'_sites.'+ctl[0].split('.')[1]        
        treefile_line = ctl[1].split('= ')[0]+'= '+str(gene)+'tree'+ctl[1].split('tree')[-1]
        ctl_new = [seqfile_line, treefile_line]
        ctl_new += ctl[2:]
        
        outf = open('baseml.ctl', 'w')
        for line in ctl_new:
            outf.write(line)
        outf.close()

def traverse_tree(gene, t_ids_gene, sf): # sf - file with sites
    '''
    Opens 'rst'+str(gene)+'sites', traverses tree, puts names&site into the nodes. Then the same with rst and genes, returns the tree with sites and needed synonymous positions from genes.
    '''

    # for the sites    
    inf = open(sf, 'r')
    a = ''.join(inf.readlines())
    inf.close()
    seqs = a.split('List of extant and reconstructed sequences')[1].split('Overall accuracy of')[0].split('\n')[4:-3]
    # expression to get the site
    dic_seqs = {'_'.join(i.split('  ')[0].split()): ''.join([j for j in i.split('  ')[1:] if j != ''][0].strip().split()) for i in seqs}

    c = 1 # counter for internal nodes
    names = []
    for leaf in t_ids_gene.traverse("preorder"):
        if leaf.name == '':
            leaf.name = 'node_#'+str(len(t_ids_gene) + c)
            c+=1
            names.append(leaf.name)
            leaf.add_features(site = dic_seqs[leaf.name])
            continue
        if leaf.name != '':
            names.append(leaf.name)
            leaf.add_features(site = dic_seqs[leaf.name])
            
    return t_ids_gene

# to get concatenated ancestral sites
def get_concat_anc_sites(tf_dic_cons, sf):

    trial_sites = open(sf, 'r')
    all_ = ''.join(trial_sites.readlines())
    trial_sites.close()
    probs = all_.split('Prob of best state at each node, listed by site')[1].split('Summary of changes along branches.')[0].strip()
    
    list_loci = []
    for l in probs.split('\n')[2:]:
        letter_prob_list = l.split(':')[1].strip()
        states = ''
        for state in letter_prob_list.split():
            states += state.strip()[0]
        list_loci.append(states)
    cons_loci = ''.join([e for i, e in enumerate(list_loci) if tf_dic_cons[i+1] != '']) # dic_cons
    return cons_loci

def fill_nodes(t, pwm, dic_cons, inc_threshold = 0):
    '''
    add the difference in weight with a parent to each node
    and a dict with nonconsensus like this {'16A': [0, 0, 0], '12C': [0, 0, 0]...}
    first number - distance to the parent, if the parent has the same nonconsensus
    second number - distance to the parent, if the energy was preserved given that nonconsensus was preserved
    third number - continuous distnace in nodes to the node where this nonconsensus first appeared
    '''
    for node in t.traverse("preorder"):
        node.add_features(weight = round(np.sum([pwm.at[e.lower(), i] for i, e in enumerate(node.site)]), 4))

    all_nc_raw = []

    for node in t.traverse("preorder"):
        if node.name == t.name:
            node.add_features(w_increment = 0)
            nc_dic = {}
            for i,e in enumerate(node.site):
                if e != dic_cons[tf][i+1] and dic_cons[tf][i+1] != '': #!!!
                    nc_dic[str(i+1)+(e)] = [0, 0, 0]
            node.add_features(Nc_dic = nc_dic)
        else:
            node.add_features(w_increment = node.weight - node.up.weight)
            nc_dic = {}
            for i,e in enumerate(node.site):
                if e != dic_cons[tf][i+1] and dic_cons[tf][i+1] != '': #!!!
                    nc_dic[str(i+1)+(e)] = [0, 0, 0]
                    all_nc_raw.append(str(i+1)+(e))
            node.add_features(Nc_dic = nc_dic)

    all_nc = list(set(all_nc_raw))

    # fill in Nc_dic
    for node in t.traverse("preorder"):
        if node.name != t.name:
            for k in node.Nc_dic.keys():
                if k in node.up.Nc_dic.keys():
                    node.Nc_dic[k][0] = round(node.dist, 2)
                    node.Nc_dic[k][2] = node.up.Nc_dic[k][2] + 1
                    if np.abs(node.w_increment) <= inc_threshold:
                        node.Nc_dic[k][1] = round(node.dist, 2)
    return t, all_nc
            
# calculate lenrat
def get_lenrat(t, trial_file_name):
    trial_file = open(trial_file_name, 'w')
    
    l_ratios_i = []
    l_ratios_i_comments = []
    
    # calculate lenrat

    for node in t.traverse("preorder"):
        if node.name != t.name:
            for k in node.Nc_dic.keys():
                if k not in node.up.Nc_dic.keys() and node.is_leaf() == False: # nonconsensus appeared?
                    l = []
                    for kid in node.traverse("preorder"): # nonconsensus evolved?
                        if kid.name != node.name:
                            if k in kid.Nc_dic.keys() and kid.Nc_dic[k][2] == node.get_distance(node, kid, topology_only=True):
                                l.append(kid.Nc_dic[k])
                    df_l = pd.DataFrame(l)
                    if df_l.empty == False:
                        length = np.sum(df_l[0])
                        length_E = np.sum(df_l[1])
                        if length != 0:
                            r = round(length_E/length, 2)
                            l_ratios_i.append(r)
                            l_ratios_i_comments.append(node.name + '__' + k)
    trial_file.close()
    
    return t, l_ratios_i, l_ratios_i_comments


# Big cycle for all TFs to get R

dic_cons = {'FadR': {1: '', 2: 'T', 3: 'C', 4: 'T', 5: 'G', 6: 'G', 7: 'T', 8: 'C', 9: '', 10: 'G', 11: 'A', 12: 'C', 13: 'C', 14: 'A', 15: 'G', 16: 'T', 17: ''},
            'HutC': {1: '', 2: '', 3: '', 4: 'C', 5: 'T', 6: 'T', 7: 'G', 8: 'T', 9: 'A', 10: 'T', 11: 'A', 12: 'T', 13: 'A', 14: 'C', 15: 'A', 16: 'A', 17: 'G', 18: '', 19: '', 20: ''},
            'HypR': {1: 'A', 2: 'T', 3: 'T', 4: 'G', 5: 'T', 6: 'A', 7: 'T', 8: 'A', 9: 'C', 10: 'A', 11: 'A', 12: 'T'},
            'PdhR': {1: 'A', 2: 'A', 3: 'T', 4: 'T', 5: 'G', 6: 'G', 7: 'T', 8: '', 9: '', 10: '', 11: 'A', 12: 'C', 13: 'C', 14: 'A', 15: 'A', 16: 'T', 17: 'T'},
            'SO0072': {1: 'T', 2: 'G', 3: 'T', 4: 'A', 5: 'T', 6: '', 7: 'A', 8: 'A', 9: 'T', 10: '', 11: '', 12: 'A', 13: 'T', 14: 'T', 15: '', 16: 'A', 17: 'T', 18: 'A', 19: 'C', 20: 'A'}, 
            'PrpR': {1: 'A', 2: 'T', 3: 'T', 4: 'G', 5: 'T', 6: 'C', 7: 'G', 8: 'A', 9: 'C', 10: 'A', 11: 'A', 12: 'T'},
            'LexA': {1: '', 2: 'A', 3: 'C', 4: 'T', 5: 'G', 6: 'T', 7: '', 8: 'T', 9: '', 10: 'T', 11: 'A', 12: '', 13: 'A', 14: '', 15: 'A', 16: 'C', 17: 'A', 18: 'G', 19: 'T', 20: ''},
            'Fnr': {1: 'T', 2: 'T', 3: 'G', 4: 'A', 5: 'T', 6: '', 7: 'T', 8: 'A', 9: '', 10: 'A', 11: 'T', 12: 'C', 13: 'A', 14: 'A'},
            'Fur': {1: '', 2: 'A', 3: 'A', 4: 'T', 5: 'G', 6: 'A', 7: '', 8: 'A', 9: 'A', 10: 'T', 11: '', 12: 'A', 13: 'T', 14: 'T', 15: '', 16: 'T', 17: 'C', 18: 'A', 19: 'T', 20: 'T', 21: ''},
            'MetJ': {1: 'A', 2: 'G', 3: 'A', 4: 'C', 5: 'G', 6: 'T', 7: 'C', 8: 'T', 9: 'A', 10: 'G', 11: 'A', 12: 'C', 13: 'G', 14: 'T', 15: 'C', 16: 'T'},
            'GcvA': {1: 'A', 2: 'T', 3: 'T', 4: 'A', 5: '', 6: '', 7: '', 8: '', 9: '', 10: '', 11: '', 12: 'T', 13: 'A', 14: 'A', 15: 'T'},
            'NarP': {1: 'T', 2: 'A', 3: 'C', 4: 'C', 5: '', 6: 'C', 7: 'T', 8: 'T', 9: 'A', 10: 'A', 11: 'G', 12: '', 13: 'G', 14: 'G', 15: 'T', 16: 'A'},
            'ArgR': {1: '', 2: '', 3: 'T', 4: 'G', 5: 'A', 6: 'A', 7: 'T', 8: '', 9: '', 10: '', 11: '', 12: 'A', 13: 'T', 14: 'T', 15: 'C', 16: 'A', 17: '', 18: ''},
            'Crp': {1: 'A', 2: 'A', 3: '', 4: 'T', 5: 'G', 6: 'T', 7: 'G', 8: 'A', 9: 'T', 10: '', 11: 'T', 12: 'A', 13: '', 14: 'A', 15: 'T', 16: 'C', 17: 'A', 18: 'C', 19: 'A', 20: '', 21: 'T', 22: 'T'},
            'FabR': {1: 'A', 2: 'G', 3: 'C', 4: '', 5: 'T', 6: 'A', 7: 'C', 8: 'A', 9: '', 10: '', 11: 'T', 12: 'G', 13: 'T', 14: 'A', 15: '', 16: 'G', 17: 'C', 18: 'T'},
            'HexR': {1: '', 2: 'T', 3: 'G', 4: 'T', 5: 'A', 6: 'A', 7: 'T', 8: '', 9: '', 10: '', 11: 'A', 12: 'T', 13: 'T', 14: 'A', 15: 'C', 16: 'A', 17: ''},
            'PsrA': {1: 'A', 2: 'T', 3: 'T', 4: '', 5: 'A', 6: 'A', 7: 'A', 8: 'C', 9: 'A', 10: '', 11: '', 12: 'T', 13: 'G', 14: 'T', 15: 'T', 16: 'T', 17: '', 18: 'A', 19: 'A', 20: 'T'},
            'TyrR': {1: '', 2: 'T', 3: 'G', 4: 'T', 5: 'A', 6: 'A', 7: 'A', 8: '', 9: '', 10: '', 11: '', 12: '', 13: '', 14: 'T', 15: 'T', 16: 'T', 17: 'A', 18: 'C', 19: 'A', 20: ''}
           }

t_ids = Tree('shew49_treeWithGenomeIds.nwk', format=2)
ancestor = t_ids.get_common_ancestor('745411.6','867.3')
t_ids.set_outgroup(ancestor)
ancestor.detach()

l_ratios_tf_list = []
l_ratios_tf_comments_list = []

inc_threshold = 0

trials = 10
problematic_trials = 100

min_spec = 9 # minimal alignment size
tfs = ['FadR', 'Crp', 'ArgR', 'FabR', 'Fnr', 'Fur', 'GcvA', 'HexR', 'LexA', 'MetJ', 'NarP', 'PsrA', 'TyrR', 'HypR', 'PdhR', 'PrpR', 'SO0072']

for tf in tfs:
    print (tf)
    dic_tree_tf = {}
    l_ratios_tf = []
    l_ratios_tf_comments = []
    # get pwm for raw sites (with paralogues and orphanes)
    f = pd.read_excel(tf+'/output.xlsx', index_col=0)
    lst = [list(str(i)) for i in f.site]
    sites = pd.DataFrame(lst, columns =[i for i in range(len(f.site[0]))])

    # make pwm
    length = sites.shape[1]
    pwm0 = pd.DataFrame(columns=[j for j in range(length)], index=['a', 'c', 'g', 't'])

    for i in range(sites.shape[1]):
        pwm0[i] = [list(sites[i]).count('a')+1, list(sites[i]).count('c')+1, list(sites[i]).count('g')+1, list(sites[i]).count('t')+1]
    # pwm normalization
    pwm = (pwm0/(sites.shape[0]+4)).applymap(lambda x: math.log2(x))
    pwm = pwm*-1
    for i in pwm.columns:
        pwm[i] = pwm[i] - np.min(pwm[i])
    pwm.to_csv(tf+'/pwm_norm')

    # get genes list
    f = pd.read_excel(tf + '/clustered.xlsx', index_col = [0])

    min_spec = 9 # minimal alignment size
    genes = list(set([i for i in list(f.index.dropna()) if list(f.index.dropna()).count(i) > min_spec]))
    print ('genes', genes)
    
    problematic_genes = []
    
    for gene in genes:
        print ('gene ', gene)
        # files gene_prot.fasta, gene.fasta , gene_sites.fasta are already there after 4nes.py script
        # making a file .nuc
        bacs = get_nuc(tf, gene)
        get_nuc_sites(tf, gene)
        # pruning a tree just for the species with this gene
        t_ids_gene = t_ids.copy('cpickle')
        t_ids_gene.prune(bacs)
        print (t_ids_gene)
        t_ids_gene.write(format=1, outfile=str(gene)+'tree.nwk')

        # open pre-created in 4nes.py true_sites_file
        true_sites_file = f'{tf}/true_sites_file_folder/true_sites_file_{gene}'
                            
        # getting a tree with sites in the nodes
        t = traverse_tree(gene, t_ids_gene, true_sites_file)
        # fill in the tree nodes
        t, nc_all = fill_nodes(t, pwm, dic_cons, inc_threshold = 0)
        # calculate r = length_E/length
        t, l_ratios, l_ratios_comments = get_lenrat(t, 'trial_file')
        dic_tree_tf[gene] = t
        l_ratios_tf.extend(l_ratios)
        l_ratios_tf_comments.extend(l_ratios_comments)          
                    
        # deleting files
        tree_file = str(gene) + 'tree.nwk'
        bashCommand = f'rm {tree_file}'
        os.system(bashCommand)

    # summing up for the tf
    l_ratios_tf_list.append(l_ratios_tf)
    l_ratios_tf_comments_list.append(l_ratios_tf_comments)
    
    # plotting
    plt.rcParams['font.size'] = '14'
    plt.figure(figsize=(6, 4))

    plt.hist(l_ratios_tf, bins = 50, color='#228B22', alpha=.6) #, density = True)
    plt.xlabel('R')
    plt.ylabel('count')
    plt.xlim(0,1)
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.title(f'Distribution of lengths ratio, {tf}')
    plt.savefig(f'{tf}/R.png', dpi=400, transparent=True)
    plt.savefig(f'pics/{tf}_R.png', dpi=400, transparent=True)

# plot for all tfs
l_ratios_tf_list_merged = []
for l_ratios_tf in l_ratios_tf_list:
    l_ratios_tf_list_merged.extend(l_ratios_tf)

plt.rcParams['font.size'] = '14'
plt.figure(figsize=(6, 4))

plt.hist(l_ratios_tf_list_merged, bins = 50, color='#228B22', alpha=.6) #, density = True)
plt.xlabel('R')
plt.ylabel('count')
plt.xlim(0,1)
plt.style.use('seaborn-v0_8-whitegrid')
plt.title(f'Distribution of lengths ratio, Shewanellaceae')
plt.savefig('pics/R_shew.png', dpi=400, transparent=True)

