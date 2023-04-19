####Download GTF file from GENCODE 
Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf:
	curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
	gzip -d gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
	mkdir Annotation
	mv gencode.v19.chr_patch_hapl_scaff.annotation.gtf Annotation/

####Processing GTF file (I type)

Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf: Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
	grep 'protein_coding\|nonsense_mediated_decay' Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf | \
	awk '{print$$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"substr($$12,2,length($$12)-3)}' > Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf

####Processing GTF file (II type)

Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID_gene_type.gtf : Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
	grep 'protein_coding\|nonsense_mediated_decay' Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf | \
	awk '{print$$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"substr($$12,2,length($$12)-3)"\t"substr($$18,2,length($$18)-3)"\t"substr($$20,2,length($$20)-3)}' > Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID_gene_type.gtf


nmd_or_protein_coding_exons.bed : Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf
	awk '{if($$3 == "exon") print$$0}' Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf | sort -u > nmd_or_protein_coding_exons.gtf
	gff2bed < nmd_or_protein_coding_exons.gtf > nmd_or_protein_coding_exons.bed

nmd_or_protein_coding_stops.bed : Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf
	awk '{if($$3 == "stop_codon") print$$0}' Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID.gtf | sort -u > nmd_or_protein_coding_stops.gtf	
	gff2bed < nmd_or_protein_coding_stops.gtf > nmd_or_protein_coding_stops.bed	

#### Intersect all stop codons with it's own exon from the same Transcripd-id

intersect/bedtools_intersect_exons_stops_1 : nmd_or_protein_coding_stops.bed nmd_or_protein_coding_exons.bed
	mkdir intersect
	bedtools intersect -wa -wb -a nmd_or_protein_coding_exons.bed -b nmd_or_protein_coding_stops.bed > intersect/bedtools_intersect_exons_stops
	awk '{if($$10 == $$20) print $$0}' intersect/bedtools_intersect_exons_stops > intersect/bedtools_intersect_exons_stops_1
	#rm nmd_or_protein_coding_exons.bed nmd_or_protein_coding_stops.bed nmd_or_protein_coding_stops.gtf nmd_or_protein_coding_exons.gtf intersect/bedtools_intersect_exons_stops

#### Intersect only nmd stop codons with it's own exon from the same Transcript-id

intersect/bedtools_intersect_exons_stops_1_nmd : Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
	awk '{if($$3 == "stop_codon") print $$0}' Annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf | \
	grep 'transcript_type "nonsense_mediated_decay"' | awk '{print$$1"\t"$$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8"\t"substr($$12,2,length($$12)-3)}' | \
	sort -u > nmd_stop_codons.gtf 
	awk 'NR == FNR{a[$$9];next} $$10 in a' nmd_stop_codons.gtf intersect/bedtools_intersect_exons_stops_1 > intersect/bedtools_intersect_exons_stops_1_nmd
	rm nmd_stop_codons.gtf		

##Covariate counting

Scripts/covar_count/exon_lenght_sum/exon_lenght_sum-after_stop_all_transcript.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/exon_lenght_sum/exon_lenght_sum-after_stop.py

Scripts/covar_count/exon_lenght_sum/exon_lenght_sum_before_PTC_all_nmd_transcript.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/exon_lenght_sum/exon_lenght_sum_before_stop.py

Scripts/covar_count/exon_lenght_sum/intron_sum_after_PTC_all_transcript.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/intron_lenght_sum/intron_sum_after_stop.py

Scripts/covar_count/exon_lenght_sum/intron_sum_before_PTC_all_transcript.txt : intersect/bedtools_intersect_exons_stops_1	
	python3 Scripts/covar_count/intron_lenght_sum/intron_sum_before_stop.py	

Scripts/covar_count/length_to_end_exon_contain_stop/length_to_exon_end_containing_PTC_all_NMD_transcript.txt : intersect/bedtools_intersect_exons_stops_1	
	python3 Scripts/covar_count/length_to_end_exon_contain_stop/length_to_exon_end_containing_stop.py

Scripts/covar_count/length_to_first_SJ/length_to_first_SJ_all_NMD_transcript_strand_sel.txt : intersect/bedtools_intersect_exons_stops_1	
	python3 Scripts/covar_count/length_to_first_SJ/length_to_first_SJ.py

Scripts/covar_count/length_to_last_SJ/length_to_last_SJ.txt : intersect/bedtools_intersect_exons_stops_1	
	python3 Scripts/covar_count/length_to_last_SJ/length_to_last_SJ.py

Scripts/covar_count/length_to_last_SJ/length_to_exon_start_containing_PTC_all_NMD_transcript.txt : intersect/bedtools_intersect_exons_stops_1	
	python3 Scripts/covar_count/length_to_start_exon_contain_stop/length_to_exon_start_containing_stop.py

Scripts/covar_count/num_exons/exons_num_after_stop.txt : intersect/bedtools_intersect_exons_stops_1		
	python3 Scripts/covar_count/num_exons/exons_num_after_stop.py	

Scripts/covar_count/num_exons/exons_num_before_stop.txt : intersect/bedtools_intersect_exons_stops_1		
	python3 Scripts/covar_count/num_exons/exons_num_before_stop.py

Scripts/covar_count/relative_stop_position_in_exon/relative_PTC_position_in_exon.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/relative_stop_position_in_exon/relative_stop_position_in_exon.py	

Scripts/covar_count/relative_stop_position_in_gene/relative_PTC_postion_in_gene.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/relative_stop_position_in_gene/relative_stop_postion_in_gene.py

Scripts/covar_count/relative_PTC_position_in_transcript/relative_PTC_position_in_transcript.txt : intersect/bedtools_intersect_exons_stops_1
	python3 Scripts/covar_count/relative_PTC_position_in_transcript/relative_PTC_position_in_transcript.py

#### Intersect poison exon stop's with its exon

## Get poison transcript ids

Scripts/covar_count/Poison_exon/Poison_id_only.txt : intersect/bedtools_intersect_exons_stops_1 Annotation/gencode.v19.chr_patch_hapl_scaff.annotation_ID_gene_type.gtf
	python3 Scripts/covar_count/Poison_exon/Poison_exon_finder.py
	awk '{if($$2==1) print $$1}' Scripts/covar_count/Poison_exon/Poison_id.txt > Scripts/covar_count/Poison_exon/Poison_id_only.txt

## Get intersect with poison transcript id

intersect/id_transID_poison.txt: intersect/bedtools_intersect_exons_stops_1
	awk 'NR == FNR{a[$$1];next} $$10 in a' Scripts/Poison/Poison_id_only.txt intersect/bedtools_intersect_exons_stops_1 > intersect/bedtools_intersect_exons_stops_1_poison
	awk '{print $$1"_"$$2+1"_"$$3"_"$$6"\t"$$10}' intersect/bedtools_intersect_exons_stops_1_poison > intersect/id_transID_poison.txt

## Get pic

pic_final : 
	python3 Scripts/Pic/exon_lenght_sum/exon_lenght_sum-after_PTC.py
	python3 Scripts/Pic/exon_lenght_sum/exon_lenght_sum-before_PTC.py
	python3 Scripts/Pic/gc_content/gc_content_after_stop.py
	python3 Scripts/Pic/gc_content/gc_content_before_stop.py
	python3 Scripts/Pic/intron_lenght_sum/intron_lenght_after_stop.py
	python3 Scripts/Pic/intron_lenght_sum/intron_lenght_before_stop.py
	python3 Scripts/Pic/length_to_end_exon_contain_stop/length_to_exon_end_containing_stop.py
	python3 Scripts/Pic/length_to_first_SJ/lenght_to_first_SJ.py
	python3 Scripts/Pic/length_to_last_SJ/lenght_to_last_SJ.py
	python3 Scripts/Pic/length_to_start_exon_contain_stop/length_to_exon_start_containing_stop.py
	python3 Scripts/Pic/num_exons/exon_num_downstream.py
	python3 Scripts/Pic/num_exons/exon_num_upstream.py
	python3 Scripts/Pic/relative_stop_position_in_exon/relative_stop_position_in_exon.py
	python3 Scripts/Pic/relative_stop_position_in_gene/relative_stop_position_in_gene.py
	python3 Scripts/Pic/relative_stop_position_in_transcript/relative_PTC_position_in_transcript.py
	mkdir Pic
	mv *.png ./Pic

pic: intersect/bedtools_intersect_exons_stops_1 intersect/bedtools_intersect_exons_stops_1_nmd intersect/id_transID_poison.txt pic_final