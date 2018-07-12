
more gdc_sample_sheet.2018-07-13.tsv |cut -f1,2,7|perl -pe 's/\t/\//g'|perl -pe 's/gz\//gz\t/g' > sample.txt
more sample.annotation_tp53_kras_braf.txt |cut -f1|perl -pe 's/\-/\t/g'|cut -f1-4|perl -pe 's/\t/\-/g' > 450k_annotation_tcga.txt

more sample.txt |grep -f 450k_annotation_tcga.txt

more sample.annotation_tp53_kras_braf.txt |cut -f1|perl -pe 's/\-/\t/g'|cut -f1-4|perl -pe 's/\t/\-/g' | \
paste - sample.annotation_tp53_kras_braf.txt |cut -f1,4 > 450k_annotation_tcga_p53.txt
