#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/num_exons


import pandas as pd

intersect = pd.read_csv('../../intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('../../Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1
intersect[20] = (intersect[12] + intersect[11])/2
intersect[20] = intersect[20].astype('int64')

def exons_num_before_PTC(exons_df,start,end):
    PTC_index = exons_df[(exons_df['start'] == start) & (exons_df['end'] == end)].index[0]
    num_exons_after_PTC = PTC_index
    return num_exons_after_PTC


for index, i in intersect.iterrows():
    file_for_write = open('./exon_num_before_PTC_all_NMD_transcript.txt','a')
    ID_transcript = i[9]
    start_exon = i[1]
    end_exon = i[2]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    select = exons_num_before_PTC(sel_exons,start_exon,end_exon)
    file_for_write.write('{}\t{}\n'.format(ID_transcript,select))
    file_for_write.close()
