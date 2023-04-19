import pandas as pd

intersect = pd.read_csv('intersect/bedtools_intersect_exons_stops_1', sep='\t', header=None, names=[i for i in range(20)])
cols = ['chr_num', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation = pd.read_csv('Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf', sep='\t', names=cols)
annotation['lenght'] = annotation['end'] - annotation['start']
intersect[11] = intersect[11] + 1
intersect[1] = intersect[1] + 1

for index, i  in intersect.iterrows():
    file = open('length_to_exon_end_containing_PTC_all_NMD_transcript.txt','a')
    if i[5] == '+':
        len = i[2] - i[12]
    else:
        len = i[11] - i[1]
    #print('{}\t{}\n'.format())
    file.write('{}\t{}\n'.format(i[9],len))
    file.close()
