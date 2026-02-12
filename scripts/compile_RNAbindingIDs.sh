#!/usr/bin/env bash
# This is a script that compiles all RNA binding IDs (InterPro entries) into a file
# and then uses that file to subset the human protein annotations to just RNA binding domains.
# The RNA binding IDs are saved to RNAbinding_IDs.tsv and the annotation subset to protein_annotation_human_rbd.tsv.

# to extract RNA recognition motifs, sort, deduplicate and write to new file to store RNA binding data
grep "RNA rec" data/protein_annotation_human.tsv | awk -F '\t' 'BEGIN{OFS="\t"}{print $2, $3}' > data/RNAbinding_IDs.tsv

# to extract RNA[ -]binding, sort, deduplicate and append to RNA binding file
grep -E "RNA[ -]binding" data/protein_annotation_human.tsv | awk -F '\t' 'BEGIN{OFS="\t"}{print $2, $3}' >> data/RNAbinding_IDs.tsv

# to extract RRMs, sort, deduplicate and append to RNA binding file
grep "RRM" data/protein_annotation_human.tsv | awk -F '\t' 'BEGIN{OFS="\t"}{print $2, $3}' >> data/RNAbinding_IDs.tsv

# to extract DNA/RNA[ -]binding and append to RNA binding file
grep -E "DNA/RNA[ -]binding" data/protein_annotation_human.tsv | awk -F '\t' 'BEGIN{OFS="\t"}{print $2, $3}' >> data/RNAbinding_IDs.tsv

# to extract all other InterPro entries listed under RNA binding superfamily IPR035979 and append to RNA binding file
grep -E "IPR005120|IPR034447|IPR002343|IPR006548|IPR027157|IPR050374|IPR050907|IPR051974|IPR012677" data/protein_annotation_human.tsv | awk -F '\t' 'BEGIN{OFS="\t"}{print $2, $3}' >> data/RNAbinding_IDs.tsv

# to extract all RNA bindings of GO terminology and append to RNA binding file
awk -F '\t' 'BEGIN{OFS="\t"}/^InterPro:IPR[0-9]+/ && /RNA binding/ {line=substr($0,10); split(line,a," > "); split(a[1],b," "); id=b[1]; desc=substr(a[1],length(id)+2); gsub(/^[ \t]+|[ \t]+$/, "", id); gsub(/^[ \t]+|[ \t]+$/, "", desc); print id "\t" desc}' data/interpro2go >> data/RNAbinding_IDs.tsv

# to remove duplicates in RNA binding file
sort data/RNAbinding_IDs.tsv | uniq > tmp && mv tmp data/RNAbinding_IDs.tsv

# to create a human protein annotation subset of just RNA binding domains based on the RNA binding file
awk -F '\t' 'BEGIN{OFS="\t"}NR==FNR {ids[$1]; next} $2 in ids' data/RNAbinding_IDs.tsv data/protein_annotation_human.tsv > data/protein_annotation_human_rbd.tsv
