#!/usr/bin/env python
#PBS -l walltime=15:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/exon_lenght_sum

import pandas as pd

intersect = pd.read_csv('intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1


def exon_lenght_sum_after_PTC(exons_df,strand,PTC_start,PTC_end,exon_contain_PTC_start,exon_contain_PTC_end, ID_transcript):
    exon_index_contain_PTC = exons_df[(exons_df['start'] == exon_contain_PTC_start) & (exons_df['end'] == exon_contain_PTC_end)].index.values[0]
    len = 0
    if strand == '+':
        len += exon_contain_PTC_end - PTC_end
    else:
        len += PTC_start - exon_contain_PTC_start
    if exon_index_contain_PTC < exons_df.shape[0] - 1:    
        exons_after_PTC_df = exons_df[exon_index_contain_PTC + 1:]
        len += (exons_after_PTC_df['end'] - exons_after_PTC_df['start']).sum()    
    return len




for index, i in intersect.iterrows():
    file_for_write = open('./exon_lenght_sum-after_stop_all_transcript.txt','a')
    ID_transcript = i[9]
    strand = i[5]
    exon_contain_PTC_start = i[1]
    exon_contain_PTC_end = i[2]
    PTC_start = i[11]
    PTC_end = i[12]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    #print(sel_exons)
    lenght = exon_lenght_sum_after_PTC(sel_exons,strand,PTC_start,PTC_end,exon_contain_PTC_start,exon_contain_PTC_end, ID_transcript)
    file_for_write.write('{}\t{}\t{}\n'.format(ID_transcript, lenght, strand))
    file_for_write.close()
