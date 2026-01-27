# Script to get figures of patterns of selection from substitution counts (run after 4nes.py)

import pandas as pd
import numpy as np
import math
from Bio import SeqIO
import csv
from Bio.Seq import Seq
from ete3 import Tree
from Bio import AlignIO
import sys
import matplotlib
from matplotlib import mlab
import matplotlib.pyplot as plt
import seaborn as sns
import os

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

#pc_fix = 0.000001

tfs = ['Crp', 'ArgR', 'FabR', 'Fnr', 'Fur', 'GcvA', 'HexR', 'LexA', 'MetJ', 'NarP', 'PsrA', 'TyrR', 'HypR', 'PdhR', 'PrpR', 'SO0072', 'FadR']
for tf in tfs:

    print (tf)

    dt_c_n = pd.read_excel(tf+'/df_c_n.xlsx', index_col=0)
    dt_n_c = pd.read_excel(tf+'/df_n_c.xlsx', index_col=0)
    dt_c_n['position'] = dt_c_n['position'].apply(lambda x: str(int(x)))
    dt_n_c['position'] = dt_n_c['position'].apply(lambda x: str(int(x)))
    pos_order = [str(x) for x in dic_cons[tf].keys()]

    sns.set(font_scale=1.3)
    sns.set_style('white')

    # plotting for CN->NCN

    # calculate pseudocount
    n = []
    for i in dt_c_n[dt_c_n['group'] == 'expected'][dt_c_n['c->n'] == 0].groupby(['gene']):
        n.append(i[1].shape[0])
    pc = 1/np.mean(n)

#     # for fixed pseudocount pc_fix
#     pc = pc_fix

    print('pc for c->n', pc)

    def find_ratio_c_n(df):
        return df[df['c->n'] == 1].shape[0]/df['c->n'].shape[0]

    grouping = dt_c_n.groupby(['position', 'group', 'gene']).apply(find_ratio_c_n)
    g = pd.DataFrame(grouping)
    g = g.reset_index()
    g['frequency of substitution'] = g[0]
    print(g)
    g_c_n = g.drop(0, axis=1)

    # add pc to obs & exp frequencies for every gene
    g_c_n['frequency of substitution'] = g_c_n['frequency of substitution'] + np.array([pc]*g_c_n.shape[0])

    # add a row with expected to the positions without such a row
    gr = g_c_n.groupby(['position', 'group']).mean()
    grp = pd.DataFrame(gr)
    grp = grp.drop('frequency of substitution', axis=1)
    pos_list = list(set(grp.index.get_level_values(0)))
    for p in pos_list:
        if grp.loc[p].shape[0] == 1:
            gene = int(grp.loc[p]['gene']) # gene doesn't matter, just to gene number
            ind = g_c_n.shape[0]
            g_c_n.loc[ind] = {'position': p, 'group': 'expected', 'gene': gene, 'frequency of substitution': pc}
    g_c_n.to_excel(tf+'/g_c_n.xlsx')

    # plotting
    c_n = sns.catplot(x='position', y='frequency of substitution', hue='group', data=g_c_n, 
                    hue_order=['observed', 'expected'], errorbar=('pi', 95),
                    order=pos_order, height=6, kind='bar', palette={'observed': '#2e8b57', 'expected': 'skyblue'})
    c_n.despine(left=True)
    c_n.set_ylabels('frequency of substitution')
    c_n.fig.suptitle(tf+' Cn->Nc', fontsize=18)
    plt.grid()
    c_n.savefig(tf+'/c_n.png', dpi=500, transparent=True)

    # plotting for n->c

    # calculate pc
    n = []
    for i in dt_n_c[dt_n_c['group'] == 'expected'][dt_n_c['n->c'] == 0].groupby(['gene']):
        n.append(i[1].shape[0])
    pc = 1/np.mean(n)

#     # for fixed pseudocount pc_fix
#     pc = pc_fix

    print('pc for n->c', pc)

    def find_ratio_n_c(df):
        return df[df['n->c'] == 1].shape[0]/df['n->c'].shape[0]

    grouping = dt_n_c.groupby(['position', 'group', 'gene']).apply(find_ratio_n_c)
    g = pd.DataFrame(grouping)
    g = g.reset_index()
    g['frequency of substitution'] = g[0]
    g_n_c = g.drop(0, axis=1)

    # add pc to obs & exp frequencies for every gene
    g_n_c['frequency of substitution'] = g_n_c['frequency of substitution'] + np.array([pc]*g_n_c.shape[0])

    # add a row with expected to the positions without such a row
    gr = g_n_c.groupby(['position', 'group']).mean()
    grp = pd.DataFrame(gr)
    grp = grp.drop('frequency of substitution', axis=1)
    pos_list = list(set(grp.index.get_level_values(0)))
    for p in pos_list:
        if grp.loc[p].shape[0] == 1:
            gene = int(grp.loc[p]['gene']) # gene doesn't matter, just to gene number
            ind = g_n_c.shape[0]
            g_n_c.loc[ind] = {'position': p, 'group': 'expected', 'gene': gene, 'frequency of substitution': pc}
    g_n_c.to_excel(tf+'/g_n_c.xlsx')

    # plotting
    n_c = sns.catplot(x='position', y='frequency of substitution', hue='group', data=g_n_c, 
                      hue_order=['observed', 'expected'], errorbar=('pi', 95),
                      order=pos_order, height=6, kind='bar', palette={'observed': '#fa8072', 'expected': 'skyblue'})
    n_c.despine(left=True)
    n_c.set_ylabels('frequency of substitution')
    n_c.fig.suptitle(tf+' Nc->Cn', fontsize=18)
    plt.grid()
    n_c.savefig(tf+'/n_c.png', dpi=500, transparent=True)

    # for 4Nes
    dt_log = pd.DataFrame(columns=['position', 'value', 'c or n'])

    # calculating just a ratio (no ln)
    def make_4nes_hm(obs_value, position, mean_exp):
        exp_value = float(mean_exp.loc[str(position),'expected'])
        nes_value = obs_value/exp_value
        return nes_value

    # with ln
    def make_4nes(obs_value, position, mean_exp):
        exp_value = float(mean_exp.loc[str(position),'expected'])
        nes_value = math.log(obs_value/exp_value)
        return nes_value

    # for CN->NCN

    # getting mean_exp with mean exp values
    grouping = g_c_n.groupby(['position', 'group']).mean()
    mean_exp = pd.DataFrame(grouping).drop('gene', axis=1)
    # getting g_c_n just with observed values
    g_c_n_observed = g_c_n[g_c_n['group'] == 'observed']

    # adding to dt_log
    for i in range(g_c_n_observed.shape[0]):
        obs_value = float(g_c_n_observed.iloc[i]['frequency of substitution'])
        position = str(g_c_n_observed.iloc[i]['position'])
        nes_value = make_4nes(obs_value, position, mean_exp)
        ind = dt_log.shape[0]
        dt_log.loc[ind] = {'position': position, 'value': nes_value, 'c or n': 'Cn->Nc'}

    # for NC->CNC

    # getting mean_exp with mean exp values
    grouping = g_n_c.groupby(['position', 'group']).mean()
    mean_exp = pd.DataFrame(grouping).drop('gene', axis=1)
    # getting g_c_n just with observed values
    g_n_c_observed = g_n_c[g_n_c['group'] == 'observed']

    # adding to dt_log
    for i in range(g_n_c_observed.shape[0]):
        obs_value = float(g_n_c_observed.iloc[i]['frequency of substitution'])
        position = str(g_n_c_observed.iloc[i]['position'])
        nes_value = make_4nes(obs_value, position, mean_exp)
        ind = dt_log.shape[0]
        dt_log.loc[ind] = {'position': position, 'value': nes_value, 'c or n': 'Nc->Cn'}

    dt_log.to_excel(tf+'/dt_log.xlsx')

    # to get sub index
    def get_sub(x):
        normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
        sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
        res = x.maketrans("".join(normal), "".join(sub_s))
        return x.translate(res)

    # plotting
    sns.set(font_scale=3)
    sns.set_style('white')

    nes = sns.catplot(x='position', y='value', hue='c or n', data=dt_log, order=pos_order, errorbar=('pi', 95),
                    height=9, aspect=2.5, kind='bar', palette={'Cn->Nc': '#2e8b57', 'Nc->Cn': '#fa8072'})
    nes.despine(left=True)
    nes.set_ylabels('4N{}s'.format(get_sub('e')))
    nes.set(ylim=(-9, 7))
    num_pos = len(dic_cons[tf].keys())
    nes.set(xlim=(-1, num_pos+0.5))
    sns.move_legend(nes, loc='center right', title='')
    new_labels = [r'CN$\rightarrow$NCN', r'NCN$\rightarrow$CN']
    for t, l in zip(nes._legend.texts, new_labels):
        t.set_text(l)
    plt.grid()
    nes.savefig(tf+'/nes.png', dpi=400, transparent=True)
    nes.savefig('pics/'+tf+'_nes.png', dpi=400, transparent=True)

    # to get violin plots
    vio_nes = sns.catplot(x='position', y='value', hue='c or n', hue_order=['Cn->Nc','Nc->Cn'], data=dt_log, order=pos_order, kind='violin',
                    width=1, height=7, aspect=3.3, palette={'Cn->Nc': '#2e8b57', 'Nc->Cn': '#fa8072'})
    vio_nes.despine(left=True)
    vio_nes.set_ylabels('4N{}s'.format(get_sub('e')))
    nes.set(ylim=(-9, 7))
    num_pos = len(dic_cons[tf].keys())
    vio_nes.set(xlim=(-1, num_pos+0.5))
    sns.move_legend(vio_nes, loc='center right', title='')
    new_labels = [r'CN$\rightarrow$NCN', r'NCN$\rightarrow$CN']
    for t, l in zip(vio_nes._legend.texts, new_labels):
        t.set_text(l)
    plt.grid()
    vio_nes.savefig(tf+'/vio_nes.png', dpi=400, transparent=True)
    vio_nes.savefig('pics/'+tf+'_vio_nes.png', dpi=400, transparent=True)

    # to get box plots
    box_nes = sns.catplot(x='position', y='value', hue='c or n', hue_order=['Cn->Nc','Nc->Cn'], data=dt_log, order=pos_order, kind='box',
                    width=1, height=7, aspect=3.3, palette={'Cn->Nc': '#2e8b57', 'Nc->Cn': '#fa8072'})
    box_nes.despine(left=True)
    box_nes.set_ylabels('4N{}s'.format(get_sub('e')))
    nes.set(ylim=(-9, 7))
    num_pos = len(dic_cons[tf].keys())
    box_nes.set(xlim=(-1, num_pos+0.5))
    sns.move_legend(box_nes, loc='center right', title='')
    new_labels = [r'CN$\rightarrow$NCN', r'NCN$\rightarrow$CN']
    for t, l in zip(box_nes._legend.texts, new_labels):
        t.set_text(l)
    plt.grid()
    box_nes.savefig(tf+'/box_nes.png', dpi=400, transparent=True)
    box_nes.savefig('pics/'+tf+'_box_nes.png', dpi=400, transparent=True)



