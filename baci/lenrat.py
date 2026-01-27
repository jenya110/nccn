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
dic_patric0 = {'Caldalkalibacillus_thermarum_TA2.A1': '986075.14',
 'Salicibibacter_kimchii_strain_NKC1-1': '2099786.3',
 'Salicibibacter_cibi_strain_NKC21-4': '2743001.3',
 'Salicibibacter_halophilus_strain_NKC3-5': '2502791.3',
 'Salicibibacter_cibarius_strain_NKC5-3': '2743000.3',
 'Paraliobacillus_zengyii_strain_X-1125': '2213194.3',
 'Radiobacillus_deserti_strain_TKL69': '2594883.4',
 'Amphibacillus_xylanus_NBRC_15112': '698758.3',
 'Oceanobacillus_zhaokaii_strain_160': '2052660.3',
 'Psychrobacillus_glaciei_strain_PB01': '2283160.3',
 'Salimicrobium_jeotgali_strain_MJ3': '1230341.3',
 'Lentibacillus_amyloliquefaciens_strain_LAM0015': '1472767.6',
 'Anaerobacillus_isosaccharinicus_strain_NB2006': '1532552.12',
 'Anoxybacillus_gonensis_strain_G2': '198467.4',
 'Anoxybacillus_caldiproteolyticus_strain_U458': '247480.5',
 'Bacillus_vini_strain_JCM_19841': '1476025.3',
 'Anoxybacillus_amylolyticus_strain_DSM_15939': '294699.3',
 'Parageobacillus_toebii_strain_NEB718': '153151.40',
 'Geobacillus_subterraneus_strain_CPW16': '129338.18',
 'Virgibacillus_dokdonensis_strain_21D': '302167.6',
 'Weizmannia_coagulans_DSM_1_=_ATCC_7050': '1121088.10',
 'Terribacillus_goriensis_strain_MP602': '386490.3',
 'Bacillus_smithii_strain_DSM_4216': '1479.5',
 'Parageobacillus_caldoxylosilyticus_strain_ER4B': '81408.26',
 'Oceanobacillus_kimchii_X50': '1238184.3',
 'Cytobacillus_gottheilii_strain_1839': '859144.10',
 'Virgibacillus_phasianinus_strain_LM2416': '2017483.3',
 'Geobacillus_zalihae_strain_SURF-189': '1220596.5',
 'Alkalihalobacillus_miscanthi_strain_AK13': '2598861.3',
 'Caldibacillus_thermoamylovorans_strain_SSBM': '35841.11',
 'Geobacillus_thermocatenulatus_strain_KCTC_3921': '33938.10',
 'Virgibacillus_halodenitrificans_strain_PDB-F2': '1482.5',
 'Virgibacillus_necropolis_strain_LMG_19488': '163877.4',
 'Halobacillus_mangrovi_strain_KTB_131': '402384.3',
 'Bacillus_amyloliquefaciens_IT-45': '1091041.4',
 'Halobacillus_halophilus_DSM_2266': '866895.3',
 'Peribacillus_psychrosaccharolyticus_strain_FDAARGOS_1161': '1407.5',
 'Bacillus_tequilensis_strain_EA-CB0015': '227866.22',
 'Virgibacillus_pantothenticus_strain_DSM_26': '1473.18',
 'Bacillus_halotolerans_strain_ZB201702': '260554.51',
 'Lysinibacillus_fusiformis_strain_RB-21': '28031.4',
 'Lysinibacillus_sphaericus_strain_DSM_28': '1421.62',
 'Fictibacillus_phosphorivorans_strain_G25-29': '1221500.4',
 'Parageobacillus_thermoglucosidasius_strain_Wild_Type': '1426.44',
 'Lysinibacillus_parviboronicapiens_strain_VT1065': '436516.5',
 'Cytobacillus_kochii_strain_BDGP4': '859143.5',
 'Peribacillus_butanolivorans_strain_PHB-7a': '421767.6',
 'Cytobacillus_ciccensis_strain_5L6': '1670641.3',
 'Bacillus_circulans_strain_FDAARGOS_783': '1397.31',
 'Mesobacillus_jeotgali_strain_DSM_18226': '129985.4',
 'Sutcliffiella_cohnii_strain_DSM_6307': '33932.3',
 'Bacillus_luti_strain_FJ': '2026191.6',
 'Bacillus_licheniformis_strain_SCDB_14': '1402.136',
 'Bacillus_cereus_strain_BC33': '1396.2858',
 'Priestia_aryabhattai_strain_K13': '412384.30',
 'Peribacillus_simplex_NBRC_15720_=_DSM_1321': '1349754.5'}

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
    
    
    sep_lines = a.split('\n')[:-1]
    bac_names = [i[1:] for i in sep_lines[::2]] # species names
    sites = sep_lines[1::2]

    bacs = [] # numeric ids
    for l in bac_names:
        bacs.append(dic_patric0[l])
    long_bac = max([len(bac) for bac in bacs]) + 3
        
    codon_file = open(tf+'/'+str(gene)+'_sites.nuc', 'w')
    codon_file.write('  '+str(num_seq)+'   '+str(len_seq)+'\n')
    for bac, site in zip(bacs, sites):
        codon_file.write(bac+' '*(long_bac-len(bac))+site+'\n')
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
#     print(list_loci)
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
            
# to calculate lenrat    
def get_lenrat(t, trial_file_name):
    trial_file = open(trial_file_name, 'w')
    
    l_ratios_i = []
    l_ratios_i_comments = []
    
    # calculating lenrat

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


dic_cons = {'CcpA': {1: '', 2: 'T', 3: 'G', 4: '', 5: 'A', 6: 'A', 7: '', 8: 'C', 9: 'G', 10: '', 11: 'T', 12: 'T', 13: '', 14: 'C', 15: 'A', 16: ''},
            'CodY': {1: 'A', 2: 'A', 3: 'T', 4: 'T', 5: 'T', 6: 'T', 7: '', 8: '', 9: '', 10: 'A', 11: 'A', 12: 'A', 13: 'A', 14: 'T', 15: 'T'}, 
            'Fur': {1: 'A', 2: '', 3: 'T', 4: 'G', 5: 'A', 6: '', 7: 'A', 8: 'A', 9: 'T', 10: '', 11: 'A', 12: 'T', 13: 'T', 14: '', 15: 'T', 16: 'C', 17: 'A', 18: '', 19: 'T'},
            'LexA': {1: '', 2: '', 3: 'G', 4: 'A', 5: 'A', 6: 'C', 7: 'A', 8: 'T', 9: 'A', 10: 'T', 11: 'G', 12: 'T', 13: 'T', 14: 'C', 15: '', 16: ''},
            'PerR': {1: 'T', 2: 'T', 3: 'A', 4: 'T', 5: 'A', 6: 'A', 7: 'T', 8: '', 9: 'A', 10: 'T', 11: 'T', 12: 'A', 13: 'T', 14: 'A', 15: 'A'},
            'PurR': {1: 'A', 2: 'A', 3: 'A', 4: '', 5: '', 6: '', 7: 'G', 8: 'A', 9: 'A', 10: '', 11: 'A', 12: '', 13: 'T', 14: '', 15: '', 16: '', 17: '', 18: '', 19: '', 20: '', 21: '', 22: '', 23: '', 24: '', 25: '', 26: '', 27: '', 28: '', 29: '', 30: '', 31: '', 32: 'A', 33: '', 34: 'T', 35: '', 36: 'T', 37: 'T', 38: 'C', 39: '', 40: '', 41: '', 42: 'T', 43: 'T', 44: 'T'},
            'TnrA': {1: 'T', 2: 'G', 3: 'T', 4: '', 5: 'A', 6: '', 7: 'A', 8: '', 9: '', 10: '', 11: 'T', 12: '', 13: 'T', 14: '', 15: 'A', 16: 'C', 17: 'A'}}

t_ids = Tree('baci58_output_treeWithGenomeIds.nwk', format=2)
Staphylococcus = t_ids.search_nodes(name='596317.3')[0]
Listeria = t_ids.search_nodes(name='1002366.3')[0]
Staphylococcus.detach()
t_ids.set_outgroup(Listeria)
Listeria.detach()

min_spec = 9 # minimal alignment size
tfs = ['Fur', 'PerR', 'PurR', 'LexA', 'TnrA', 'CodY', 'CcpA']

l_ratios_tf_list = []
l_ratios_tf_comments_list = []

inc_threshold = 0

trials = 10
problematic_trials = 100

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

    genes = list(set([i for i in list(f.index.dropna()) if list(f.index.dropna()).count(i) > min_spec]))
    print ('genes', genes)
    
    problematic_genes = []
    
    for gene in genes:
        print ('gene ', gene)
        # files gene_prot.fasta, gene.fasta , gene_sites.fasta are already there after 4nes.py script
#        get_nt(tf, gene, f) # if not, make them again
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
    # write l_ratios_tf in the file
    l_ratios_tf_file = open(f'{tf}/l_ratios_tf.txt', 'w')
    l_ratios_tf_file.write(str(l_ratios_tf))
    l_ratios_tf_file.close()
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

# write l_ratios for all tfs in the file
l_ratios_tf_file = open('l_ratios_all_tf.txt', 'w')
l_ratios_tf_file.write(str(l_ratios_tf_list_merged))
l_ratios_tf_file.close()

plt.rcParams['font.size'] = '14'
plt.figure(figsize=(6, 4))

plt.hist(l_ratios_tf_list_merged, bins = 50, color='#228B22', alpha=.6) #, density = True)
plt.xlabel('R')
plt.ylabel('count')
plt.xlim(0,1)
plt.style.use('seaborn-v0_8-whitegrid')
plt.title(f'Distribution of lengths ratio, Shewanellaceae')
plt.savefig('pics/R_shew.png', dpi=400, transparent=True)

