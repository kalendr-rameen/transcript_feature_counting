#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/intron_lenght_sum

import pandas as pd


intersect = pd.read_csv('intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1


def intron_sum_before_PTC(exons_df,exon_contain_PTC_start,exon_contain_PTC_end,strand):
    exon_contain_PTC_index = exons_df[(exons_df['start'] == exon_contain_PTC_start) & (exons_df['end'] == exon_contain_PTC_end)].index[0]
    if strand == '+':
        exons_len_sum = exons_df.iloc[:exon_contain_PTC_index]['lenght'].sum()
        intron_sum = (exon_contain_PTC_start - exons_df.iloc[[0]]['start'].values[0]) - exons_len_sum
    else:
        exons_len_sum = exons_df.iloc[:exon_contain_PTC_index]['lenght'].sum()
        intron_sum = (exons_df.iloc[[0]]['end'].values[0] - exon_contain_PTC_end) - exons_len_sum
        
    return intron_sum


#file_for_write = open('./dists5.txt','a') 
for index, i in intersect.iterrows():
    file_for_write = open('./intron_sum_before_PTC_all_transcript.txt','a')
    ID_transcript = i[9]
    exon_contain_PTC_start = i[1]
    exon_contain_PTC_end = i[2]
    strand = i[5]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    select = intron_sum_before_PTC(sel_exons,exon_contain_PTC_start,exon_contain_PTC_end,strand)
    file_for_write.write('{}\t{}\t{}\n'.format(ID_transcript,select,strand))
    file_for_write.close()
#file.close()