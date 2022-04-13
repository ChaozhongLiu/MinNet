for i in $(seq 1 8)
do
    Rscript liger.R $i
    python3 after_liger.py $i
done
