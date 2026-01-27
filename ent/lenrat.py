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
dic_patric0 = {'Citrobacter_rodentium_NBRC_105723_=_DSM_16636_chromosome': '67825.26',
 'Klebsiella_sp._PCX_chromosome': '2686364.3',
 'Raoultella_sp._HC6_chromosome': '2923366.4',
 'Buttiauxella_ferragutiae_strain_H4-C11_chromosome': '82989.7',
 'Kluyvera_ascorbata_strain_SK_chromosome': '51288.34',
 'Lelliottia_sp._AC1_chromosome': '2067959.4',
 'Cedecea_lapagei_strain_NCTC11466_chromosome_1': '158823.5',
 'Citrobacter_youngae_strain_NCTC13709_chromosome_1': '133448.8',
 'Klebsiella_grimontii_strain_4928STDY7071328_chromosome_1': '571.706',
 'Pseudescherichia_vulneris_strain_4928STDY7071515_chromosome_1': '566.7',
 'Enterobacter_sp._BIDMC_29_strain_Colony302_chromosome': '1329841.4',
 'Citrobacter_sp._Colony219_chromosome': '2861802.3',
 'Escherichia_sp._TM-G17TGC_chromosome': '2860337.3',
 'Enterobacter_cancerogenus_strain_JY65_chromosome': '69218.53',
 'Klebsiella_sp._CTHL.F3a_chromosome': '2873296.5',
 'Kluyvera_sp._CRP_chromosome': '2873269.5',
 'Enterobacter_sp._JBIWA003_chromosome': '2831890.5',
 'Salmonella_sp._A29-2_chromosome': '2876620.6',
 'Klebsiella_variicola_subsp._variicola_strain_F2R9T_chromosome': '2590157.4',
 'Klebsiella_variicola_subsp._tropicalensis_strain_CDC4241-71': '2489014.3',
 'Klebsiella_quasivariicola_strain_08A119_chromosome': '2026240.29',
 'Plesiomonas_shigelloides_strain_7A_chromosome_1': '703.94',
 'Leclercia_sp._G3L_chromosome': '2898725.3',
 'Pseudocitrobacter_sp._G163CM_chromosome': '2891570.3',
 'Klebsiella_sp._KPN54798_chromosome': '2897413.6',
 'Klebsiella_pasteurii_strain_Kox205_chromosome': '2587529.4',
 'Citrobacter_sp._MGH_55_strain_Colony470_chromosome': '1439319.4',
 'Citrobacter_sp._BIDMC107_strain_Colony252_chromosome': '1686384.4',
 'Leclercia_pneumoniae_strain_49125_chromosome': '2815358.3',
 'Salmonella_sp._SAL-007_chromosome': '2816952.3',
 'Enterobacter_sp._BWH_37_chromosome': '1329836.3',
 'Citrobacter_sedlakii_strain_3347689II_chromosome': '67826.48',
 'Klebsiella_sp._P1927_chromosome': '2829836.3',
 'Klebsiella_sp._A52_chromosome': '2834819.3',
 'Citrobacter_pasteurii_strain_FDAARGOS_1424_chromosome': '1563222.16',
 'Citrobacter_farmeri_strain_FDAARGOS_1423_chromosome': '67824.35',
 'Salmonella_enterica_subsp._arizonae_strain_LHICA_AZ23_chromosome': '59203.506',
 'Salmonella_enterica_subsp._salamae_strain_LHICA_SA1_chromosome': '59202.401',
 'Lelliottia_amnigena_strain_FDAARGOS_1445_chromosome': '61646.100',
 'Escherichia_sp._TC-EC600-tetX4_chromosome': '2857061.3',
 'Kosakonia_sp._SMBL-WEM22_chromosome': '2725560.3',
 'Blochmannia_endosymbiont_of_Colobopsis_nipponica_isolate_CNIPBac': '2681987.3',
 'Enterobacter_bugandensis_strain_STN0717-56_chromosome': '881260.71',
 'Klebsiella_michiganensis_strain_THO-011_chromosome': '1134687.230',
 'Kosakonia_pseudosacchari_strain_BDA62-3_chromosome': '1646340.7',
 'Citrobacter_sp._BDA59-3_chromosome': '2781952.3',
 'Klebsiella_sp._BDA134-6_chromosome': '2787706.3',
 'Enterobacter_mori_strain_HSW1412_chromosome': '539813.12',
 'Enterobacter_sp._YSU_chromosome': '648691.3',
 'Klebsiella_quasipneumoniae_strain_KqPF26_chromosome': '1463165.411',
 'Klebsiella_africana_strain_FF1003_chromosome': '2489010.3',
 'Citrobacter_sp._R56_chromosome': '1573676.3',
 'Citrobacter_cronae_strain_Colony478_chromosome': '1748967.19',
 'Escherichia_albertii_strain_Sample_167_chromosome': '208962.170',
 'Cedecea_sp._FDAARGOS_727_chromosome': '2545798.3',
 'Escherichia_sp._SCLE84_chromosome': '2725997.3',
 'Phytobacter_diazotrophicus_strain_UAEU22_chromosome': '395631.4',
 'Kosakonia_sp._MUSA4_chromosome': '2067958.3',
 'Salmonella_sp._SCFS4_chromosome': '2725417.3',
 'Salmonella_enterica_subsp._VII_serovar_1': '41520.3',
 'Buttiauxella_agrestis_strain_DSM_9389_chromosome': '82977.7',
 'Shigella_sonnei_strain_SE6-1_chromosome': '624.2249',
 'Citrobacter_sp._172116965_chromosome': '2683822.3',
 'Klebsiella_sp._RHBSTW-00464_chromosome': '2742649.3',
 'Citrobacter_sp._RHBSTW-00017_chromosome': '2742629.3',
 'Citrobacter_sp._RHB25-C09_chromosome': '2742624.3',
 'Citrobacter_sp._RHB20-C15_chromosome': '2742620.3',
 'Escherichia_fergusonii_strain_RHB19-C05_chromosome': '564.119',
 'Klebsiella_sp._WP3-S18-ESBL-05_chromosome': '2675711.3',
 'Kluyvera_genomosp._3_strain_PO2S7_chromosome': '2705459.3',
 'Raoultella_terrigena_strain_JH01_chromosome': '577.26',
 'Enterobacter_asburiae_strain_1808-013_chromosome': '61645.239',
 'Enterobacter_sp._18A13_chromosome': '2565914.3',
 'Klebsiella_aerogenes_strain_Ka37751_chromosome': '548.605',
 'Atlantibacter_hermannii_strain_ATCC_33651_chromosome': '565.13',
 'Kosakonia_radicincitans_strain_DSM_107547_chromosome': '283686.18',
 'Citrobacter_portucalensis_strain_FDAARGOS_617_chromosome': '1639133.19',
 'Citrobacter_werkmanii_strain_FDAARGOS_616_chromosome': '67827.25',
 'Cronobacter_sp._JZ38_chromosome': '1906275.8',
 'Enterobacter_sichuanensis_strain_SGAir0282_chromosome': '550.2510',
 'Enterobacter_oligotrophicus_strain_CCA6_chromosome': '2478464.3',
 'Kosakonia_arachidis_strain_KACC_18508_chromosome': '551989.5',
 'Citrobacter_braakii_strain_MiY-A_chromosome': '57706.65',
 'Klebsiella_variicola_strain_LEMB11_chromosome': '244366.435',
 'Citrobacter_sp._H12-3-2_chromosome': '2664853.3',
 'Kluyvera_intermedia_strain_N2-1_chromosome': '61648.20',
 'Leclercia_sp._J807_chromosome': '2681307.3',
 'Leclercia_sp._119287_chromosome': '2681308.3',
 'Leclercia_sp._1106151_chromosome': '2681309.4',
 'Citrobacter_sp._LUTT5_chromosome': '2697370.3',
 'Citrobacter_sp._Y3_chromosome': '2716879.3',
 'Enterobacter_huaxiensis_strain_090008_=_WCHEHu090008_chromosome': '2339234.3',
 'Buttiauxella_sp._3AFRM03_chromosome': '2479367.3',
 'Citrobacter_freundii_strain_FDAARGOS_549_chromosome': '546.470',
 'Klebsiella_sp._FDAARGOS_511_chromosome': '2488567.3',
 'Klebsiella_oxytoca_strain_FDAARGOS_500_chromosome': '571.376',
 'Kosakonia_sp._CCTCC_M2018092_chromosome': '2492396.3',
 'Metakosakonia_sp._MRY16-398_chromosome': '2487150.3',
 'Scandinavium_goeteborgense_strain_CCUG_66741_chromosome': '1851514.3',
 'Salmonella_sp._SSDFZ54_chromosome': '2500542.3',
 'Klebsiella_sp._LY_chromosome': '2015795.3',
 'Kosakonia_cowanii_strain_FBS_223_chromosome': '208223.19',
 'Citrobacter_sp._ABFQG_chromosome': '2529121.3',
 'Citrobacter_arsenatis_strain_LY-1_chromosome': '2546350.3',
 'Citrobacter_tructae_strain_SNU_WT2_chromosome': '2562449.3',
 'Citrobacter_sp._TBCP-5362_chromosome': '2576406.3',
 'Jejubacter_calystegiae_strain_KSNA2_chromosome': '2579935.4',
 'Escherichia_sp._E4742_chromosome': '2044467.5',
 'Citrobacter_sp._CF971_chromosome': '2566012.3',
 'Raoultella_electrica_strain_DSM_102253_chromosome': '1259973.4',
 'Klebsiella_sp._LTGPAF-6F_chromosome': '1905288.3',
 'Enterobacter_chengduensis_strain_WCHECl-C4_=_WCHECh050004': '550.1095',
 'Klebsiella_sp._M5al_chromosome': '1934254.3',
 'Klebsiella_quasipneumoniae_subsp._similipneumoniae_strain_G747': '1463164.76',
 'Cedecea_neteri_strain_FDAARGOS_392_chromosome': '158822.21',
 'Shigella_dysenteriae_strain_BU53M1_chromosome': '622.81',
 'Klebsiella_pneumoniae_subsp._ozaenae_strain_WCHKP030925_chromosome': '574.12',
 'Enterobacter_sp._Crenshaw_chromosome': '1977566.3',
 'Escherichia_marmotae_strain_HT073016_chromosome': '1499973.6',
 'Pluralibacter_gergoviae_strain_FDAARGOS_186_chromosome': '61647.65',
 'Lelliottia_sp._WB101_chromosome': '2153385.3',
 'Phytobacter_sp._SCO41_chromosome': '1756993.3',
 'Citrobacter_sp._CRE-46_strain_AR_0157_chromosome': '1703250.3',
 'Klebsiella_huaxiensis_strain_WCHKl090001_chromosome': '2153354.3',
 'Klebsiella_quasipneumoniae_subsp._quasipneumoniae_strain_A708': '1667327.19',
 'Raoultella_sp._X13_chromosome': '2259647.3',
 'Leclercia_sp._W6_chromosome': '2282310.3',
 'Cronobacter_sakazakii_strain_CS-931_chromosome': '28141.658',
 'Shigella_boydii_strain_ATCC_9210_chromosome': '621.85',
 'Citrobacter_sp._MGH103_chromosome': '1686378.3',
 'Cronobacter_universalis_NCTC_9529_chromosome': '1074000.5',
 'Cronobacter_muytjensii_ATCC_51329_chromosome': '1159613.6',
 'Cronobacter_malonaticus_LMG_23826_chromosome': '1159491.6',
 'Cronobacter_dublinensis_subsp._dublinensis_LMG_23823_chromosome': '1159554.7',
 'Cronobacter_condimenti_1330_strain_LMG_26250_chromosome': '1073999.7',
 '[Enterobacter]_lignolyticus_strain_G5_chromosome': '1715259.3',
 'Leclercia_adecarboxylata_strain_USDA-ARS-USMARC-60222_chromosome': '83655.10',
 'Citrobacter_amalonaticus_strain_FDAARGOS_165_chromosome': '35703.29',
 'Citrobacter_sp._FDAARGOS_156_strain_FDAARGOS_155_chromosome': '1702170.3',
 'Salmonella_enterica_subsp._diarizonae_strain_11-01855_chromosome': '59204.12',
 'Kosakonia_oryzae_strain_Ola_51_chromosome': '497725.5',
 'Kosakonia_sacchari_strain_BO-1_chromosome': '1158459.5',
 'Enterobacter_hormaechei_subsp._oharae_strain_DSM_16687_chromosome': '301102.37',
 'Enterobacter_roggenkampii_strain_DSM_16690_chromosome': '1812935.7',
 'Enterobacter_ludwigii_strain_EN-119_chromosome': '299767.18',
 'Plautia_stali_symbiont_DNA': '1560356.5',
 'Shigella_sp._PAMC_28760_chromosome': '1813821.3',
 'Lelliottia_jeotgali_strain_PFL01': '1907578.3',
 'Citrobacter_telavivensis_strain_6105_chromosome': '2653932.3',
 'Citrobacter_sp._S39_chromosome': '2660638.3',
 'Salmonella_sp._HNK130_chromosome': '2664291.3',
 'Salmonella_sp._S13_chromosome': '2686305.3',
 'Citrobacter_europaeus_strain_FDAARGOS_1490_chromosome': '1914243.25',
 'Salmonella_sp._JXY0409-18_chromosome': '2879113.3',
 'Klebsiella_pneumoniae_strain_12_chromosome': '573.24247',
 'Salmonella_enterica_strain_NCTC9948_genome_assembly': '28901.4088',
 'Raoultella_ornithinolytica_strain_NCTC8846_genome_assembly': '548.410',
 'Escherichia_coli_str._K-12_substr._MG1655': '511145.183',
 'Shigella_flexneri_2a_str._301_chromosome': '198214.7',
 'Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._LT2': '90371.4767',
 'Citrobacter_koseri_ATCC_BAA-895': '290338.8',
 'Enterobacter_hormaechei': '718254.4',
 'Enterobacter_soli': '640513.3',
 'Klebsiella_pneumoniae_subsp._pneumoniae_HS11286_chromosome': '1125630.4',
 'Shimwellia_blattae_DSM_4481_=_NBRC_105725': '630626.3',
 'Enterobacter_kobei': '1211025.3',
 'Salmonella_bongori_N268-08': '1197719.3',
 'Enterobacter_hormaechei_subsp._hoffmannii_ECNIH3_chromosome': '1333851.3',
 'Enterobacter_cloacae_strain_GGT036_chromosome': '550.150',
 'Raoultella_planticola_strain_FDAARGOS_64_chromosome': '575.6',
 'Enterobacter_hormaechei_subsp._xiangfangensis_strain_34978': '550.253',
 'Enterobacter_hormaechei_subsp._steigerwaltii_strain_34998': '550.254',
 'Blochmannia_endosymbiont_of_Polyrhachis_turneri_strain': '1505596.4',
 'Blochmannia_endosymbiont_of_Camponotus_obliquus_strain': '1505597.4',
 'Phytobacter_ursingii_strain_CAV1151_chromosome': '61648.5',
 'Pasteurella_multocida_P2100': '747.1971',
 'Vibrio_cholerae_N1252': '666.6203'}

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
            print('was blank')
            leaf.name = 'node_#'+str(len(t_ids_gene) + c)
            c+=1
            print ('leaf.name', leaf.name)
            names.append(leaf.name)
            leaf.add_features(site = dic_seqs[leaf.name])
            print ('leaf.site', leaf.site)
            continue
        if leaf.name != '':
            print('was not blank')
            print ('leaf.name', leaf.name)
            names.append(leaf.name)
            leaf.add_features(site = dic_seqs[leaf.name])
            print ('leaf.site', leaf.site)
            
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
    and a dictionary with nonconsensus like this {'16A': [0, 0, 0], '12C': [0, 0, 0]...}
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
            
# to calculate lenght ratio R   
def get_lenrat(t, trial_file_name): # 

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

dic_cons = {'ArgR': {1: '', 2: '', 3: 'T', 4: 'G', 5: '', 6: 'A', 7: 'T', 8: '', 9: '', 10: '', 11: '', 12: 'A', 13: 'T', 14: '', 15: 'C', 16: 'A', 17: '', 18: ''},
            'Crp': {1: '', 2: '', 3: '', 4: 'T', 5: 'G', 6: 'T', 7: 'G', 8: 'A', 9: 'T', 10: '', 11: '', 12: '', 13: '', 14: 'A', 15: 'T', 16: 'C', 17: 'A', 18: 'C', 19: 'A', 20: '', 21: '', 22: ''},
            'FadR': {1: 'A', 2: '', 3: 'C', 4: '', 5: 'G', 6: 'G', 7: 'T', 8: 'C', 9: '', 10: 'G', 11: 'A', 12: 'C', 13: 'C', 14: '', 15: 'G', 16: '', 17: 'T'},
            'Fnr': {1: 'T', 2: 'T', 3: 'G', 4: 'A', 5: 'T', 6: '', 7: 'T', 8: 'A', 9: '', 10: 'A', 11: 'T', 12: 'C', 13: 'A', 14: 'A'},
            'FruR': {1: '', 2: '', 3: 'T', 4: 'G', 5: 'A', 6: 'A', 7: '', 8: 'C', 9: 'G', 10: '', 11: 'T', 12: 'T', 13: 'C', 14: 'A', 15: '', 16: ''},            
            'Fur': {1: '', 2: '', 3: '', 4: 'A', 5: 'A', 6: 'T', 7: '', 8: 'A', 9: 'T', 10: '', 11: 'A', 12: 'T', 13: '', 14: 'A', 15: 'T', 16: 'T', 17: '', 18: '', 19: ''},
            'KdgR': {1: '', 2: '', 3: '', 4: 'T', 5: '', 6: 'A', 7: 'A', 8: 'A', 9: '', 10: '', 11: '', 12: '', 13: '', 14: 'T', 15: 'T', 16: 'T', 17: '', 18: 'A', 19: '', 20: '', 21: ''},            
            'LexA': {1: '', 2: 'A', 3: 'C', 4: 'T', 5: 'G', 6: 'T', 7: '', 8: 'T', 9: '', 10: '', 11: '', 12: '', 13: 'A', 14: '', 15: 'A', 16: 'C', 17: 'A', 18: 'G', 19: 'T', 20: ''},            
            'MetJ': {1: 'A', 2: 'G', 3: 'A', 4: '', 5: '', 6: 'T', 7: '', 8: 'T', 9: 'A', 10: '', 11: 'A', 12: '', 13: '', 14: 'T', 15: 'C', 16: 'T'},            
            'MetR': {1: '', 2: 'T', 3: 'G', 4: 'A', 5: 'A', 6: '', 7: '', 8: '', 9: '', 10: '', 11: 'T', 12: 'T', 13: 'C', 14: 'A', 15: ''},            
            'NarP': {1: 'T', 2: 'A', 3: 'C', 4: '', 5: '', 6: '', 7: '', 8: '', 9: '', 10: '', 11: '', 12: '', 13: '', 14: 'G', 15: 'T', 16: 'A'},            
            'NtrC': {1: '', 2: 'G', 3: 'C', 4: 'A', 5: 'C', 6: '', 7: '', 8: '', 9: '', 10: '', 11: '', 12: '', 13: 'G', 14: 'T', 15: 'G', 16: 'C', 17: ''},            
            'PurR': {1: 'A', 2: '', 3: 'G', 4: '', 5: 'A', 6: 'A', 7: '', 8: 'C', 9: 'G', 10: '', 11: 'T', 12: 'T', 13: '', 14: 'C', 15: '', 16: 'T'},
            'TyrR': {1: '', 2: 'T', 3: 'G', 4: 'T', 5: 'A', 6: 'A', 7: 'A', 8: '', 9: '', 10: '', 11: '', 12: '', 13: '', 14: 'T', 15: 'T', 16: 'T', 17: 'A', 18: 'C', 19: 'A', 20: ''}}

t_ids = Tree('ent_cut_tree.nwk', format=2)
Pasteurella = t_ids.search_nodes(name='747.1971')[0]
Vibrio = t_ids.search_nodes(name='666.6203')[0]
ancestor = t_ids.get_common_ancestor(Pasteurella, Vibrio)
t_ids.set_outgroup(ancestor)
ancestor.detach()

l_ratios_tf_list = []
l_ratios_tf_comments_list = []

inc_threshold = 0

trials = 10
problematic_trials = 100

min_spec = 9 # minimal alignment size
tfs = ['ArgR', 'FadR', 'Fnr', 'FruR', 'Fur', 'KdgR', 'LexA', 'MetJ', 'MetR', 'NarP', 'NtrC', 'PurR', 'TyrR', 'Crp']

for tf in tfs:
    print (tf)
    dic_tree_tf = {}
    l_ratios_tf = []
    l_ratios_tf_comments = []
    # get pwm for raw sites (with paralogues and orphanes)
    f = pd.read_excel(tf+'/output_cut.xlsx', index_col=0)
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
    f = pd.read_excel(tf + '/clustered_cut.xlsx', index_col = [0])

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
plt.title(f'Distribution of lengths ratio, Enterobacteriaceae')
plt.savefig('pics/R_ent.png', dpi=400, transparent=True)
