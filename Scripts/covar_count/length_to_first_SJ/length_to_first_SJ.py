#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/length_to_first_SJ

import pandas as pd

intersect = pd.read_csv('../../intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('../../Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1
intersect[20] = ((intersect[11] + intersect[12])/2).astype('int32')


def lenght_to_first_SJ(exons_df,strand,PTC_start,PTC_end,exon_contain_PTC_start,exon_contain_PTC_end, ID_transcript):
    exon_index_contain_PTC = exons_df[(exons_df['start'] == exon_contain_PTC_start) & (exons_df['end'] == exon_contain_PTC_end)].index.values[0]
    len = 0
    if exon_index_contain_PTC > 0:
        if strand == '+':
            len += PTC_start - exon_contain_PTC_start
        else:
            len += exon_contain_PTC_end - PTC_end

        exons_before_PTC_df = exons_df[1:exon_index_contain_PTC]
        len += (exons_before_PTC_df['end'] - exons_before_PTC_df['start']).sum()
    return len



for index, i in intersect.iterrows():
    file_for_write = open('./length_to_first_SJ_all_NMD_transcript_strand_sel.txt','a')
    ID_transcript = i[9]
    strand = i[5]
    PTC_start=i[11]
    PTC_end=i[12]
    exon_contain_PTC_start = i[1]
    exon_contain_PTC_end = i[2]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    if strand == '+':
        sel_exons.sort_values(by=['start'],ascending=True, inplace=True)
    else:
        sel_exons.sort_values(by=['start'],ascending=True, inplace=False)
    sel_exons.reset_index(inplace=True,drop=True)
    lenght = lenght_to_first_SJ(sel_exons,strand,PTC_start,PTC_end,exon_contain_PTC_start,exon_contain_PTC_end, ID_transcript)
    file_for_write.write('{}\t{}\n'.format(ID_transcript, lenght))
    file_for_write.close()
