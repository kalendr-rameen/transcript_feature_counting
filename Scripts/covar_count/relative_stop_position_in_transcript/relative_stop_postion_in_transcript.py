#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/relative_PTC_position_in_transcript

import pandas as pd

intersect = pd.read_csv('../../intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('../../Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1
intersect[20] = ((intersect[11] + intersect[12])/2).astype('int32')

def relative_position_in_transcript(exons_df,PTC_mid,strand, exon_start, exon_stop, PTC_start, PTC_stop):
    PTC_index = exons_df[(exons_df['start'] == exon_start) & (exons_df['end'] == exon_stop)].index[0]
    if strand == '+':
        value = exons_df['lenght'][:PTC_index].sum() + (PTC_mid - exon_start)
        relative_position = value / exons_df['lenght'].sum()
    else:
        value = exons_df['lenght'][:PTC_index].sum() + (exon_stop - PTC_mid)
        relative_position = value / exons_df['lenght'].sum()
    return relative_position



for index, i in intersect.iterrows():
    file_for_write = open('./relative_PTC_position_in_transcript_all_NMD_transcript.txt','a')
    ID_transcript = i[9]
    #print(ID_transcript)
    PTC_mid = i[20]
    strand = i[5]
    exon_start = i[1]
    exon_stop = i[2]
    PTC_start = i[11]
    PTC_stop = i[12]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    select = relative_position_in_transcript(sel_exons,PTC_mid,strand, exon_start, exon_stop, PTC_start, PTC_stop)
    file_for_write.write('{}\t{}\t{}\n'.format(ID_transcript,select,strand))
    file_for_write.close()
