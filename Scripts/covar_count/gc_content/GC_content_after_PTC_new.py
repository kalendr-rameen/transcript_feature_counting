#!/usr/bin/env python
#PBS -l walltime=15:00:00,mem=20gb,nodes=4:ppn=1
#PBS -d /home/rameen/Covariates/gc_content



import pandas as pd
import subprocess
import getpass
import os

intersect = pd.read_csv('/Users/ramen/Desktop/Diplom/Summer_work/Intersect/intersect_all_exons_n_stops/bedtools_intersect_exons_stops_1', sep='\t', header=None,names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'transcript_name']
annotation = pd.read_csv('/Users/ramen/Desktop/Diplom/Summer_work/GTF_file/gencode.v19.annotation-2.gtf_withproteinids_ID', sep='\t', names=cols)
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1
annotation['lenght'] = annotation['end'] - annotation['start']

intersect[20] = (intersect[12] + intersect[11])/2
intersect[20] = intersect[20].astype('int64')

def GC_content_after_PTC(exons_df,PTC_mid,strand,exon_contain_PTC_start,exon_contain_PTC_end,chr_num):
    string = ''
    file_for_write = open('example.bed', mode='w')
    exon_index_contain_PTC = exons_df[(exons_df['start'] == exon_contain_PTC_start) & (exons_df['end'] == exon_contain_PTC_end)].index.values[0]

    after_PTC_exons = exons_df.iloc[exon_index_contain_PTC + 1:]
    for index, i in after_PTC_exons.iterrows():
        file_for_write.write('{}\t{}\t{}\n'.format(i['chr_num'],i['start'],i['end']))
    if strand == '+':
        file_for_write.write('{}\t{}\t{}'.format(chr_num, PTC_mid, exon_contain_PTC_end))
    else:
        file_for_write.write('{}\t{}\t{}'.format(chr_num, exon_contain_PTC_start, PTC_mid))

    file_for_write.close()

    get_fasta = subprocess.check_output("bedtools getfasta -fi /Users/ramen/Desktop/Diplom/Summer_work/Genome/merge_genome.fasta -bed ./example.bed", shell=True)
    get_fasta_split = str(get_fasta).split('\\n')

    for i in range(1,len(get_fasta_split) - 1,2):
        string += get_fasta_split[i]
    return string


#file_for_write = open('./GC_content_after_PTC_all_nmd_transcripts.txt','a')
for index, i in intersect.iloc[73122:].iterrows():
    file_for_write = open('./GC_content_after_PTC_all_nmd_transcripts.txt','a')
    ID_transcript = i[9]
    PTC_mid = i[20]
    strand = i[5]
    exon_contain_PTC_start = i[1]
    exon_contain_PTC_end = i[2]
    chr_num = i[0]
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    select = GC_content_after_PTC(sel_exons,PTC_mid,strand,exon_contain_PTC_start,exon_contain_PTC_end, chr_num)
    try:
        file_for_write.write('{}\t{}\n'.format(ID_transcript,(select.count('G') + select.count('C'))/len(select)))
    except:
        file_for_write.write('{}\t{}\n'.format(ID_transcript,0))
    file_for_write.close()
#file_for_write.close()
