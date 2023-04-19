import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns

relative_PTC_position_in_gene_table = pd.read_csv('Scripts/covar_count/relative_stop_position_in_gene/relative_PTC_position_in_gene_all_NMD_transcript.txt', sep='\t', header=None, names=['Transcript_ID','relative_PTC_position_in_gene','strand'])
intersect = pd.read_csv('intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None,names=[i for i in range(20)])
intersect2 = pd.read_csv('intersect/bedtools_intersect_exons_stops_1_nmd', sep='\t', header=None,names=[i for i in range(20)])
intersect = pd.merge(intersect,intersect2, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
intersect.reset_index(inplace=True, drop=True)
intersect[11] = intersect[11] + 1
intersect[1] =  intersect[1] + 1
intersect.dropna(inplace=True)
new_column_1 = []
new_column_2 = []
for index, row in intersect.iterrows():
    new_column_1.append('{}_{}_{}_{}'.format(row[0],row[1],row[2],row[5]))
    new_column_2.append(row[9])
table_protein_coding = pd.DataFrame({'id': new_column_1, 'Transcript_ID': new_column_2})
table_poison = pd.read_csv('intersect/id_transID_poison.txt', sep='\t', header=None, names=['id','Transcript_ID'])
frames = [table_protein_coding, table_poison]
table = pd.concat(frames)
nmdnar_upf1 = pd.read_csv('nmdnar_data/upf1xrn1vscontrol.tsv', sep='\t')
nmdnar_smg6 = pd.read_csv('nmdnar_data/smg6xrn1vscontrol.tsv', sep='\t')
table_upf1 = pd.merge(nmdnar_upf1,table, on='id').merge(relative_PTC_position_in_gene_table, on='Transcript_ID')
table_upf1['KD'] = 'UPF1/XRN1'
table_upf1['relative_PTC_position_in_gene'] = table_upf1['relative_PTC_position_in_gene'].astype('float64')
table_smg6 = pd.merge(nmdnar_smg6,table, on='id').merge(relative_PTC_position_in_gene_table, on='Transcript_ID')
table_smg6['KD'] = 'SMG6/XRN1'
table_smg6['relative_PTC_position_in_gene'] = table_smg6['relative_PTC_position_in_gene'].astype('float64')

table_unite = pd.concat([table_smg6,table_upf1])
left_int = pd.Interval(left=-1, right=-0.1)
middle_int = pd.Interval(left=-0.1, right=0.1, closed='neither')
right_int = pd.Interval(left=0.1, right=1, closed='left')
table_unite['psi_bin_col'] = table_unite['deltaPSIc'].apply(lambda x: left_int if x <= -0.1 else (middle_int if x <= 0.1 else right_int))
table_unite.sort_values(by='psi_bin_col',inplace=True)
table_unite.rename({'psi_bin_col':'deltaPSI_бин', 'relative_PTC_position_in_gene':'относительная_позиция_стоп_кодона_в_гене', 'KD':'нокаут'}, axis=1, inplace=True)

def annotate(data, **kws):
    #print(data)
    n = len(data)
    ax = plt.gca()
    interval_1 = data['deltaPSI_бин'].unique()[0]
    interval_2 = data['deltaPSI_бин'].unique()[1]
    interval_3 = data['deltaPSI_бин'].unique()[2]
    v_1 = data[data['deltaPSI_бин'] == interval_1]['относительная_позиция_стоп_кодона_в_гене'].values
    v_2 = data[data['deltaPSI_бин'] == interval_2]['относительная_позиция_стоп_кодона_в_гене'].values
    v_3 = data[data['deltaPSI_бин'] == interval_3]['относительная_позиция_стоп_кодона_в_гене'].values
    res_1 = mannwhitneyu(v_1, v_3, method="asymptotic", alternative='two-sided')
    #res_2 = mannwhitneyu(v_2, v_3, method="asymptotic")
    #res_3 = mannwhitneyu(v_1, v_3, method="asymptotic")
    # ax.text(.9, .9, f"N = {n}", transform=ax.transAxes, fontweight='bold')
    ax.text(.9,.9, f"p-value={'{0:.4g}'.format(res_1.pvalue)}", transform=ax.transAxes, fontweight='bold')
    #ax.text(.8,.7, f"p-value={res_2.pvalue}", transform=ax.transAxes, fontweight='bold')
    #ax.text(.8,.6, f"p-value={res_3.pvalue}", transform=ax.transAxes, fontweight='bold')

s = sns.catplot(
    data=table_unite, x='относительная_позиция_стоп_кодона_в_гене', y='deltaPSI_бин',
    col='нокаут', kind='box', col_wrap=2, showfliers = False, orient='h', palette=['blue','grey','yellow']
)
s.map_dataframe(annotate)
s.set_xticklabels(rotation=30)
s.savefig('./относительная_позиция_стоп_кодона_в_гене.png')