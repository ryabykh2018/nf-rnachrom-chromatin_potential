#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate chromatin potential and generate plots')
    parser.add_argument('--input_path', type=str, required=True,
                      help='Input directory path')
    parser.add_argument('--output_path', type=str, required=True,
                      help='Output directory path')
    parser.add_argument('--input_path_RNAseq', type=str, required=True,
                      help='Input RNA-seq directory path')
    parser.add_argument('--counts_RNAseq', type=str, required=True,
                      help='RNA-seq counts filename')
    parser.add_argument('--counts_contacts', type=str, required=True,
                      help='Contacts counts filename')
    parser.add_argument('--N_contacts_min', type=int, default=100,
                      help='Minimum number of contacts (default: 100)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                      help='FDR threshold (default: 0.05)')
    parser.add_argument('--gene_len_min', type=int, default=100,
                      help='Minimum gene length (default: 100)')
    parser.add_argument('--type', type=str, required=True,
                      help='Type of contacts (UU_all or UU_dist or UM_all or UM_dist)')
    return parser.parse_args()

def calculate_stats(row, Ne, Nc):
    #P_rna_RNAseq = row['N_counts_RNAseq'] / Ne
    #P_rna_contacts = row['N_contacts'] / Nc
    #P = (row['N_counts_RNAseq'] + row['N_contacts']) / (Ne + Nc)
    #Z = (P_rna_contacts - P_rna_RNAseq) / np.sqrt(P * (1 - P) * (1/Nc + 1/Ne))
    
    # Proportion test
    Z, pval = proportions_ztest([row['N_contacts'], row['N_counts_RNAseq']], [Nc, Ne])
    return pd.Series({'chP': Z,'pval': pval})

def joinplot_chP_Ncontacts(test_data, save_name):

    test_data_sorted = test_data.sort_values('biotype2', ascending=True)
    
    edge_colors = {}
    edge_colors2 = {}
    default_colors = sns.color_palette()[:2]
    
    biotype_list = test_data_sorted['biotype'].unique()
    biotype2_dict = {}
    for i in biotype_list:
        biotype2_dict[i.split(" ")[1]] = i
        if i.split(" ")[1] == "mRNAs":
            edge_colors[i] = to_rgba(default_colors[0], 0.3)  # Синий с alpha=0.5
            edge_colors2[i] = to_rgba(default_colors[0], 1)  # Синий с alpha=0.5
        else:
            edge_colors[i] = to_rgba(default_colors[1], 0.3)   # Красный с alpha=0.5
            edge_colors2[i] = to_rgba(default_colors[1], 1)   # Красный с alpha=0.5
    
    # Создаем пустой jointplot
    g = sns.jointplot(data=test_data_sorted,x="N_contacts", y="chP",hue='biotype', alpha=0, s=50,linewidth=0,legend=False,height=7,space=0.2)
    
    # Убираем числовые метки у осей плотности
    g.ax_marg_x.set_yticklabels([])
    g.ax_marg_x.set_ylabel('')
    g.ax_marg_y.set_xticklabels([])
    g.ax_marg_y.set_xlabel('')
    # Убираем числовые значения у шкалы плотности (если она есть)
    if hasattr(g, 'ax_joint'):
        for cbar in g.fig.axes:
            if isinstance(cbar, plt.Axes) and len(cbar.get_images()) > 0:
                cbar.set_visible(False)
    
    # Применяем log-шкалу для оси X если нужно
    g.ax_joint.set_xscale('log')
    # Для логарифмической шкалы устанавливаем разумные пределы
    x_min = test_data_sorted["N_contacts"][test_data_sorted["N_contacts"] > 0].min() * 0.9
    x_max = test_data_sorted["N_contacts"].max() * 1.1
    g.ax_joint.set_xlim(x_min, x_max)

    # Ручная отрисовка точек с прозрачными контурами
    for biotype, group in test_data_sorted.groupby('biotype2'):
        g.ax_joint.scatter(
            x=group["N_contacts"],
            y=group["chP"],
            facecolor='none',          # Прозрачная заливка
            edgecolor=edge_colors[biotype2_dict[biotype]],  # Цвет контура с прозрачностью
            s=30,                      # Размер точек
            linewidth=1.8,             # Толщина контура
            label=biotype
        )
        
    # Находим точки с максимальным и предмаксимальным значениями по OY
    sorted_by_oy = test_data_sorted.sort_values("chP", ascending=False)
    max_oy_point = sorted_by_oy.iloc[0]
    second_max_oy_point = sorted_by_oy.iloc[1]
    
    # Функция для определения положения подписи
    def get_annotation_params(x_val):
        threshold = 10000
        if x_val > threshold:  # Если точка в правой части графика
            return (-10, 0), 'right'  # Подпись слева
        else:  # Если точка в левой части
            return (10, 0), 'left'  # Подпись справа
    
    # # Подпись для точки с максимальным OY 
    # xytext1, ha1 = get_annotation_params(max_oy_point['N_contacts'])
    # g.ax_joint.annotate(max_oy_point['gene_name'], xy=(max_oy_point['N_contacts'], max_oy_point["chP"]),
    #                    xytext=xytext1,textcoords='offset points',fontsize=12,ha=ha1,va='center')
    
    # # Подпись для точки с предмаксимальным OY
    # xytext2, ha2 = get_annotation_params(second_max_oy_point['N_contacts'])
    # g.ax_joint.annotate(second_max_oy_point['gene_name'], xy=(second_max_oy_point["N_contacts"], second_max_oy_point["chP"]),
    #                    xytext=xytext2,textcoords='offset points',fontsize=12,ha=ha2,va='center')

    # Настройки легенды (показываем реальные цвета)
    # legend_elements = [
    #     plt.Line2D([0], [0], marker='o', color='w', label=biotype2_dict['mRNAs'], markerfacecolor='none',markersize=10, markeredgecolor=edge_colors2[biotype2_dict['mRNAs']], markeredgewidth=1.8),
    #     plt.Line2D([0], [0], marker='o', color='w', label=biotype2_dict['ncRNAs'], markerfacecolor='none',markersize=10, markeredgecolor=edge_colors2[biotype2_dict['ncRNAs']],markeredgewidth=1.8)
    # ]
    legend_elements = []
    for biotype, label in biotype2_dict.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor='none',markersize=10, markeredgecolor=edge_colors2[label], markeredgewidth=1.8))


    g.ax_joint.legend(handles=legend_elements, title='Biotype',frameon=True, framealpha=0.5, fontsize=14, title_fontsize=16, loc='lower right')
    
    # Подписи осей
    xlabel = 'The total number of contacts'
    g.ax_joint.set_xlabel(xlabel, fontsize=16, labelpad=10)
    g.ax_joint.set_ylabel("Chromatin potential", fontsize=16, labelpad=10)
    
    # Установка размера шрифта для меток осей
    g.ax_joint.tick_params(axis='x', labelsize=14)
    g.ax_joint.tick_params(axis='y', labelsize=14)
    
    # Добавляем легкую сетку
    # g.ax_joint.grid(True, linestyle='--', alpha=0.5)
    
    # PDF (векторный формат)
    # plt.savefig(f"{save_name}.pdf", format='pdf', dpi=360, bbox_inches="tight")
    # JPG (растровый формат)
    plt.savefig(f"{save_name}.png", format='png', dpi=360, bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()
    
    # Read input files
    counts_RNAseq = pd.read_csv(os.path.join(args.input_path_RNAseq, args.counts_RNAseq), sep='\t')
    counts_RNAseq['gene_length'] = counts_RNAseq['end'] - counts_RNAseq['start'] + 1
    counts_RNAseq = counts_RNAseq[['gene_name', 'gene_type', 'N_counts', 'gene_length']].rename(columns={'N_counts': 'N_counts_RNAseq'})

    counts_contacts = pd.read_csv(os.path.join(args.input_path, args.counts_contacts), sep='\t')
    counts_contacts = counts_contacts[['gene_name', 'gene_type', 'N_counts']].rename(columns={'N_counts': 'N_contacts'})
    
    # Merge and filter data
    counts_merge = counts_RNAseq.merge(counts_contacts, how='left')
    counts_merge = counts_merge[(counts_merge['N_contacts'] > args.N_contacts_min) & (counts_merge['gene_length'] > args.gene_len_min)]
    counts_merge[['N_counts_RNAseq', 'N_contacts']] += 1
    
    if len(counts_merge) != 0:
        # Calculate statistics
        Nc = counts_merge['N_contacts'].sum()
        Ne = counts_merge['N_counts_RNAseq'].sum()
        
        chP = counts_merge.join(counts_merge.apply(lambda x: calculate_stats(x, Ne, Nc), axis=1))
        chP['fdr_bh'] = multipletests(chP['pval'], method='fdr_bh')[1]

        # Save results
        chP_output = chP.copy()
        chP_output['gene_length'] = chP_output['gene_length'].astype(int)
        chP_output['N_contacts'] = chP_output['N_contacts'].astype(int)
        chP_output['N_counts_RNAseq'] = chP_output['N_counts_RNAseq'].astype(int)
        chP_output['chP'] = chP_output['chP'].round(2)
        chP_output['pval'] = chP_output['pval'].map('{:.2e}'.format)
        chP_output['fdr_bh'] = chP_output['fdr_bh'].map('{:.2e}'.format)
        chP_output.rename(columns={'N_contacts': 'N_contacts_plus_1', 'N_counts_RNAseq': 'N_counts_RNAseq_plus_1'}).to_csv(os.path.join(args.output_path, 'chP_{type}.tab'.format(type=args.type)), sep='\t', index=False)

        if chP[chP['fdr_bh'] < args.fdr_threshold].shape[0] != 0:
            chP = chP[chP['fdr_bh'] < args.fdr_threshold].copy()
            chP.loc[chP['gene_type'] == 'protein_coding', 'biotype'] = str(chP[chP['gene_type'] == 'protein_coding'].shape[0]) + ' mRNAs'
            chP.loc[chP['gene_type'] != 'protein_coding', 'biotype'] = str(chP[chP['gene_type'] != 'protein_coding'].shape[0]) + ' ncRNAs'
            chP.loc[chP['gene_type'] == 'protein_coding', 'biotype2'] = 'mRNAs'
            chP.loc[chP['gene_type'] != 'protein_coding', 'biotype2'] = 'ncRNAs'
            joinplot_chP_Ncontacts(chP, save_name=os.path.join(args.output_path, 'chP_{type}'.format(type=args.type)))


if __name__ == '__main__':
    main()
