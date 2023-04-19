#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/relative_PTC_position_in_gene


import pandas as pd

intersect = pd.read_csv('../../intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('../../Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1
intersect[20] = ((intersect[11] + intersect[12])/2).astype('int32')

def relative_position_in_gene(exons_df,PTC_mid,strand):
   if strand == '+':
    start_pos = exons_df['start'].min()
    end_pos = exons_df['end'].max()
    relative = (PTC_mid - start_pos) / (end_pos - start_pos)
   else:
    start_pos = exons_df['start'].min()
    end_pos = exons_df['end'].max()
    relative = 1 - ((PTC_mid - start_pos) / (end_pos - start_pos))
   return relative


for index, i in intersect.iterrows():
    file_for_write = open('./relative_PTC_position_in_gene_all_NMD_transcript.txt','a')
    ID_transcript = i[9]
    #print(ID_transcript)
    PTC_mid = i[20]
    strand = i[5]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    select = relative_position_in_gene(sel_exons,PTC_mid,strand)
    file_for_write.write('{}\t{}\t{}\n'.format(ID_transcript,select,strand))
    file_for_write.close()
