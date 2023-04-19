import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import getpass
import os
import subprocess
import numpy

intersect = pd.read_csv('intersect/bedtools_intersect_exons_stops_1', sep='\t',header=None,names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'transcript_name', 'gene_name', 'transcript_type']
annotation = pd.read_csv('Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID_gene_type.gtf', sep='\t', names=cols)
temporary = annotation[['transcript_name','gene_name']]
temporary[9] = temporary['transcript_name']
temporary.drop('transcript_name',axis=1, inplace=True)
temporary.drop_duplicates(inplace=True)
intersect[9] = intersect[9].astype(object)
temporary[9] = temporary[9].astype(object)
intersect = pd.merge(temporary, intersect, on=9, how='right')
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1

def select_flank_exon_coord(sel_exons,exon_contain_PTC_start,exon_contain_PTC_end, strand):
    PTC_exon_index = sel_exons[(sel_exons['start'] == exon_contain_PTC_start) & (sel_exons['end'] == exon_contain_PTC_end)].index[0]
    
    if sel_exons.shape[0] == PTC_exon_index + 1:
        flank_left = sel_exons.iloc[PTC_exon_index - 1]
        flank_right = pd.Series(index=annotation.columns, data=[numpy.nan for x in range(sel_exons.shape[1])])
    elif PTC_exon_index == 0:
        flank_left = pd.Series(index=annotation.columns, data=[numpy.nan for x in range(sel_exons.shape[1])])
        flank_right = sel_exons.iloc[PTC_exon_index + 1]
    else:
        flank_left = sel_exons.iloc[PTC_exon_index - 1]
        flank_right = sel_exons.iloc[PTC_exon_index + 1]
    flank_coords_left = numpy.array([flank_left['start'],flank_left['end']])
    flank_coords_right = numpy.array([flank_right['start'],flank_right['end']])
    #print(len(flank_coords_left),len(flank_coords_right))
    return flank_coords_left,flank_coords_right

#flank_right = sel_exons.iloc[PTC_exon_index + 1]

def poison_exon_finder(sel_exons, coords_left, coords_right, exon_contain_PTC_start, exon_contain_PTC_end, trans_ID):
    
    transcript_ID = sel_exons['transcript_name'].drop_duplicates().values
    if numpy.isnan(coords_left).any():    # select only     
        for ID in transcript_ID:
            sel_exons_by_ID = sel_exons[sel_exons['transcript_name'] == ID]
#             sel_row_by_ID = sel_exons_by_ID[((((sel_exons_by_ID['start'] == coords_right[0])& (sel_exons_by_ID['end'] == coords_right[1])) \
#                             |  ((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end)))\
#                             & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            sel_right_exon = sel_exons_by_ID[(((sel_exons_by_ID['start'] == coords_right[0]) & (sel_exons_by_ID['end'] == coords_right[1]))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            
            sel_PTC_contain_exon = sel_exons_by_ID[((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding')]
            if sel_right_exon.shape[0] == 1 and sel_PTC_contain_exon.shape[0] == 0:
                return 1
                break
                #print(1,ID, trans_ID)
            else:
                continue
                #print(0,ID, trans_ID)
                
    if numpy.isnan(coords_right).any(): 
        for ID in transcript_ID:
            sel_exons_by_ID = sel_exons[sel_exons['transcript_name'] == ID]
#             sel_row_by_ID = sel_exons_by_ID[((((sel_exons_by_ID['start'] == coords_right[0])& (sel_exons_by_ID['end'] == coords_right[1])) \
#                             |  ((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end)))\
#                             & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            sel_left_exon = sel_exons_by_ID[(((sel_exons_by_ID['start'] == coords_left[0]) & (sel_exons_by_ID['end'] == coords_left[1]))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            
            sel_PTC_contain_exon = sel_exons_by_ID[((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding')]
            if sel_left_exon.shape[0] == 1 and sel_PTC_contain_exon.shape[0] == 0:
                return 1
                break
                #print(1,ID, trans_ID)
            else:
                continue
                #print(0,ID, trans_ID)
    else:
        for ID in transcript_ID:
            sel_exons_by_ID = sel_exons[sel_exons['transcript_name'] == ID]
#             sel_row_by_ID = sel_exons_by_ID[((((sel_exons_by_ID['start'] == coords_right[0])& (sel_exons_by_ID['end'] == coords_right[1])) \
#                             |  ((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end)))\
#                             & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            sel_right_exon = sel_exons_by_ID[(((sel_exons_by_ID['start'] == coords_right[0]) & (sel_exons_by_ID['end'] == coords_right[1]))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            sel_left_exon = sel_exons_by_ID[(((sel_exons_by_ID['start'] == coords_left[0]) & (sel_exons_by_ID['end'] == coords_left[1]))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding'))]
            sel_PTC_contain_exon = sel_exons_by_ID[((sel_exons_by_ID['start'] == exon_contain_PTC_start) & (sel_exons_by_ID['end'] == exon_contain_PTC_end))                                            & (sel_exons_by_ID['transcript_type'] == 'protein_coding')]
            if sel_left_exon.shape[0] == 1 and sel_PTC_contain_exon.shape[0] == 0 and sel_right_exon.shape[0] == 1:
                return 1
                break
                #print(1,ID, trans_ID)
            else:
                continue
                #print(0,ID, trans_ID)
    return 0

for index, i in intersect.iterrows():
    file_for_write = open('./Poison_id.txt','a')
    ID_transcript = i[9]
    #PTC_mid = i[20]
    strand = i[5]
    exon_contain_PTC_start = i[1]
    exon_contain_PTC_end = i[2]
    chr_num = i[0]
    gene_name = i['gene_name']
    sel_exons = annotation[(annotation[annotation.columns[8]] == ID_transcript) & (annotation[annotation.columns[2]] == 'exon')]
    sel_exons.reset_index(inplace=True,drop=True)
    coords_left, coords_right = select_flank_exon_coord(sel_exons,exon_contain_PTC_start,exon_contain_PTC_end, strand)
    all_exons_by_gene = annotation[(annotation['gene_name'] == gene_name) & (annotation['type'] == 'exon')]
    transcript_ID_by_gene = all_exons_by_gene['transcript_name'].drop_duplicates().values
    poison_result = poison_exon_finder(all_exons_by_gene, coords_left, coords_right, exon_contain_PTC_start, exon_contain_PTC_end, ID_transcript)                                
    file_for_write.write('{}\t{}\n'.format(ID_transcript,poison_result))