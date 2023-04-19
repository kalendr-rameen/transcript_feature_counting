#!/usr/bin/env python
#PBS -l walltime=20:00:00,mem=20gb,nodes=1:ppn=1
#PBS -d /home/rameen/Covariates/relative_PTC_position_in_exon

import pandas as pd

intersect = pd.read_csv('../../intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('../../Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1

plus_table = intersect[intersect[5] == '+']
minus_table = intersect[intersect[5] == '-']

plus_table = plus_table.assign(relative=lambda x: (x[2] - x[12])/(x[2]-x[1]))
minus_table = minus_table.assign(relative=lambda x: (x[11] - x[1])/(x[2]-x[1]))

frames = [plus_table, minus_table]
df = pd.concat(frames)
df = df.rename(columns={9: "ID_transcript", "relative": "relative_PTC_position_in_exon"})
df[['ID_transcript', 'relative_PTC_position_in_exon']].to_csv('./relative_PTC_position_in_exon_all_nmd_transcripts.csv', sep='\t', index=False)
