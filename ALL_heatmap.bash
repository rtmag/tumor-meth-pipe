
echo > /root/TCGA/TCGA/logs/tp53_HEATMAP_20-apr-2018.log

for i in ../preprocessing_tp53/*/rnb.set.norm.RData.zip;
do pathy=${i//\/rnb.set.norm.RData.zip/};
cancer=${pathy//\.\.\/preprocessing\_tp53\//};
echo $pathy;
cd $pathy ;
  Rscript /root/TCGA/TCGA/scripts/all_heatmap.R >> /root/TCGA/TCGA/logs/tp53_HEATMAP_20-apr-2018.log ; 
  new_name="/root/TCGA/TCGA/matrix/$cancer.meth.norm.sig.rds"
  mv meth.norm.sig.rds $new_name
cd ../../scripts/ ;
done


Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/ACC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/CESC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/CHOL/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/COAD/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/ESCA/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/GBM/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/HNSC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/KICH/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/KIRC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/KIRP/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/LAML/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/LGG/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/LIHC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/LUAD/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/LUSC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/MESO/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/PAAD/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/PCPG/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/PRAD/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/READ/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/SARC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/SKCM/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/STAD/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/TGCT/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/THCA/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/THYM/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/UCEC/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/UCS/ &
Rscript /root/TCGA/TCGA/scripts/all_heatmap.R  /root/TCGA/TCGA/preprocessing_tp53/UVM/ &
