
echo > /root/TCGA/TCGA/logs/tp53_HEATMAP_20-apr-2018.log

for i in ../preprocessing_tp53/*/rnb.set.norm.RData.zip;
do pathy=${i//\/rnb.set.norm.RData.zip/};
cancer=${pathy//\.\.\/preprocessing\_tp53\//};
echo $pathy;
cd $pathy ;
  Rscript /root/TCGA/TCGA/scripts/ >> /root/TCGA/TCGA/logs/tp53_HEATMAP_20-apr-2018.log ; 
cd ../../scripts/ ;
done


/root/TCGA/TCGA/matrix/
