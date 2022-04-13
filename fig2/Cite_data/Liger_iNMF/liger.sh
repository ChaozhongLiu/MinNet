for i in $(seq 1 8)
do
    Rscript Liger_iNMF.R $i
    python3 after_liger.py $i
done
