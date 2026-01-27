# Script to calculate patterns of selection from the number substitutions NCN -> CN and CN -> NCN. Files generated in this script will be needed for all the other scripts.

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
            
    # for the genes
    
    syn_codons = ['A', 'G', 'P', 'T', 'V']

    inf = open('rst', 'r')
    a = ''.join(inf.readlines())
    inf.close()
    seqs = a.split('List of extant and reconstructed sequences')[1].split('Overall accuracy of')[0].split('\n')[4:-3] # for marginal reconstruction rather than joint reconstruction
    dic_seqs = {'_'.join(i.split('  ')[0].split()): ''.join([j for j in i.split('  ')[1:] if j != ''][0].strip().split()) for i in seqs}
    
    aa = []
    nt = []
    for name in names:
        nucseq = dic_seqs[name][2::3] # only 3d nucleotides
        protseq = Seq(dic_seqs[name]).translate()
        aa.append(list(protseq))
        nt.append(list(nucseq))
    df_aa = pd.DataFrame(aa)
    df_nt = pd.DataFrame(nt)
    names_aa = [i+'_aa' for i in names]
    df_aa = df_aa.set_index([pd.Index(names_aa)])
    df_nt = df_nt.set_index([pd.Index(names)])
    df = pd.concat([df_aa,df_nt]).T
    num_nodes = 2*len(t_ids_gene)-1
    aa_patterns = [i*num_nodes for i in syn_codons]
    df['aa'] = df[names_aa].apply(lambda x: ''.join(x), axis = 1)
    df = df[df['aa'].isin(aa_patterns)]
    df = df.drop(columns=names_aa+['aa'])
    df = df.T

    for leaf in t_ids_gene.traverse("preorder"):
        leaf.add_features(gene = ''.join(list(df.loc[leaf.name])))
    return t_ids_gene


# count and get dataframe with substitutions
def count(tf, gene, t, dt_c_n, dt_n_c, dic_cons): # t - tree, dt - dataframe with substitutions used to add more rows
    for node in t.traverse("preorder"):
        if node.children != []:
            name0 = node.name
            site_seq0 = node.site
            gene_seq0 = node.gene
            
            for ch in [0, 1]: # 1st&2nd child
                name1 = node.children[ch].name
                site_seq1 = node.children[ch].site
                gene_seq1 = node.children[ch].gene
                    
                # count c->n for sites (consensus to non-consensus)
                for n in dic_cons[tf].keys():
                    if site_seq0[n-1] == dic_cons[tf][n]:
                        if site_seq0[n-1] == site_seq1[n-1]:
                            ind = dt_c_n.shape[0]
                            dt_c_n.loc[ind] = {'position': n, 'c->n': 0, 'group': 'observed', 'gene': gene, 'branch': name0+'->'+name1}
                        if site_seq0[n-1] != site_seq1[n-1]:
                            ind = dt_c_n.shape[0]
                            dt_c_n.loc[ind] = {'position': n, 'c->n': 1, 'group': 'observed', 'gene': gene, 'branch': name0+'->'+name1}
                        # count c->n for genes (consensus to non-consensus)
                        for l in range(len(gene_seq0)):
                            if gene_seq0[l] == dic_cons[tf][n]:
                                if gene_seq0[l] == gene_seq1[l]:
                                    ind = dt_c_n.shape[0]
                                    dt_c_n.loc[ind] = {'position': n, 'c->n': 0, 'group': 'expected', 'gene': gene, 'branch': name0+'->'+name1}
                                if gene_seq0[l] != gene_seq1[l]:
                                    ind = dt_c_n.shape[0]
                                    dt_c_n.loc[ind] = {'position': n, 'c->n': 1, 'group': 'expected', 'gene': gene, 'branch': name0+'->'+name1}
                # count n->c for sites (non-consensus to consensus)
                for n in dic_cons[tf].keys():
                    if site_seq0[n-1] != dic_cons[tf][n] and dic_cons[tf][n] != '':
                        if site_seq1[n-1] != dic_cons[tf][n]:
                            ind = dt_n_c.shape[0]
                            dt_n_c.loc[ind] = {'position': n, 'n->c': 0, 'group': 'observed', 'gene': gene, 'branch': name0+'->'+name1}
                        if site_seq1[n-1] == dic_cons[tf][n]:
                            ind = dt_n_c.shape[0]
                            dt_n_c.loc[ind] = {'position': n, 'n->c': 1, 'group': 'observed', 'gene': gene, 'branch': name0+'->'+name1}
                        # count n->c for genes (non-consensus to consensus)
                        for l in range(len(gene_seq0)):
                            if gene_seq0[l] == site_seq0[n-1]:
                                if gene_seq1[l] != dic_cons[tf][n]:
                                    ind = dt_n_c.shape[0]
                                    dt_n_c.loc[ind] = {'position': n, 'n->c': 0, 'group': 'expected', 'gene': gene, 'branch': name0+'->'+name1}
                                if gene_seq1[l] == dic_cons[tf][n]:
                                    ind = dt_n_c.shape[0]
                                    dt_n_c.loc[ind] = {'position': n, 'n->c': 1, 'group': 'expected', 'gene': gene, 'branch': name0+'->'+name1}
    
    return dt_c_n, dt_n_c

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


# Big cycle for all TFs to get ...df_n_c.xlsx tables

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

min_spec = 9 # minimal number of of alignments considered
tfs = ['Fur', 'PerR', 'PurR', 'LexA', 'TnrA', 'CodY', 'CcpA']

# will be needed to identify problematic genes (that give ambiguous outputs in baseml)
trials = 10
problematic_trials = 100

for tf in tfs:
    print ('tf', tf)
    
    f = pd.read_excel(tf + '/clustered.xlsx', index_col = [0])
    genes = list(set([i for i in list(f.index.dropna()) if list(f.index.dropna()).count(i) > min_spec]))
    f = f[f.index.isin(genes)]
    print (genes)

    # to make one picture for all the TF's genes
    d_c_n = {'position': [], 'c->n': [], 'group': [], 'gene': [], 'branch': []}
    dt_c_n = pd.DataFrame(data=d_c_n)
    d_n_c = {'position': [], 'n->c': [], 'group': [], 'gene': [], 'branch': []}
    dt_n_c = pd.DataFrame(data=d_n_c)

    for gene in genes:
        print ('gene ', gene)
        # making files gene_prot.fasta, gene.fasta , gene_sites.fasta
        get_nt(tf, gene, f)
        # making a file .nuc
        bacs = get_nuc(tf, gene)
        get_nuc_sites(tf, gene)
        # pruning a tree just for the species with this gene
        t_ids_gene = t_ids.copy('cpickle')
        t_ids_gene.prune(bacs)
        print (t_ids_gene)
        t_ids_gene.write(format=1, outfile=str(gene)+'tree.nwk')
        # rewriting baseml.ctl
        rewrite_ctl(tf, gene, True) # False if gene, True if site
        # running baseml for sites
        bashCommand = 'baseml'
        os.system(bashCommand)
        # moving everything from rst to в 'rst'+str(gene)+'sites'
        sites_file = f'{tf}/trials/rst'+str(gene)+'sites_'+'0'
        bashCommand = f'cp rst {sites_file}'
        os.system(bashCommand)
        
        prev_anc_sites = get_concat_anc_sites(dic_cons[tf], sites_file)
        problematic_gene = False
        true_sites_file = sites_file
        
        # to identify problematic genes (that give ambiguous outputs in baseml)
        for trial in range(1, trials):
            # running baseml for sites
            bashCommand = 'baseml'
            os.system(bashCommand)
            # moving everything from rst to в 'rst'+str(gene)+'sites'
            sites_file = f'{tf}/trials/rst'+str(gene)+'sites_'+str(trial)
            bashCommand = f'cp rst {sites_file}'
            os.system(bashCommand)
            current_anc_sites = get_concat_anc_sites(dic_cons[tf], sites_file)
            if current_anc_sites != prev_anc_sites:
                problematic_gene = True
                true_sites_file = ''
                break
        print('trial', trial)
        print('problematic gene:', problematic_gene)
        
        # to get the dominating result (winner)
        if problematic_gene == True:
            ptrials_dic = {}
            for ptrial in range(problematic_trials):
                # running baseml for sites
                bashCommand = 'baseml'
                os.system(bashCommand)
                # moving everything from rst to 'rst'+str(gene)+'sites'
                sites_file = f'{tf}/trials/rst'+str(gene)+'sites_p'+str(ptrial)
                bashCommand = f'cp rst {sites_file}'
                os.system(bashCommand)
                ptrials_dic[ptrial] = get_concat_anc_sites(dic_cons[tf], sites_file)
            concat_anc_sites = pd.Series(ptrials_dic.values())
            print(pd.DataFrame(concat_anc_sites).groupby(0).size())
            winner = concat_anc_sites.describe().top
            print(concat_anc_sites.describe())
            print('winner', winner)
            for ptrial in list(ptrials_dic.keys()):
                if ptrials_dic[ptrial] == winner:
                    true_sites_file = f'{tf}/trials/rst'+str(gene)+'sites_p'+str(ptrial)
                    break

        # rewriting codeml.ctl
        rewrite_ctl(tf, gene, False) # False if gene, True if site
        # running codeml for genes
        bashCommand = 'codeml'
        os.system(bashCommand)
        # getting a tree with sites and cut genes in nodes
        t = traverse_tree(gene, t_ids_gene, true_sites_file)
        # save the true_sites_file -- will be needed in further scripts
        bashCommand = f'cp {true_sites_file} {tf}/true_sites_file_folder/true_sites_file_{gene}'
        os.system(bashCommand)
        # moving everything from rst to 'rst'+str(gene)+'genes'
        genes_file = f'{tf}/rst'+str(gene)+'genes'
        bashCommand = f'cp rst {genes_file}'
        os.system(bashCommand)
        # deleting the file with the tree
        tree_file = str(gene) + 'tree.nwk'
        bashCommand = f'rm {tree_file}'
        os.system(bashCommand)
        # adding rows to dt_cn, dt_nc datasets
        print ('calculating...')
        dt_c_n, dt_n_c = count (tf, gene, t, dt_c_n, dt_n_c, dic_cons)
        # deleting sites files in trials/
        if problematic_gene == True:
            for ptrial in list(ptrials_dic.keys()):
                sites_file_to_delete = f'{tf}/trials/rst'+str(gene)+'sites_p'+str(ptrial)
                bashCommand = f'rm {sites_file_to_delete}'
                os.system(bashCommand)
        else:
            for t in range(trials):
                sites_file_to_delete = f'{tf}/trials/rst'+str(gene)+'sites_'+str(t)
                bashCommand = f'rm {sites_file_to_delete}'
                os.system(bashCommand)

    dt_n_c.to_excel(tf + '/df_n_c.xlsx')
    dt_c_n.to_excel(tf + '/df_c_n.xlsx')
