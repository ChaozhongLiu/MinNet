for i in $(seq 1 8)
do
    Rscript Liger.R $i
    python3 after_liger.py $i
done
