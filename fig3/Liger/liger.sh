for i in $(seq 1 8)
do
    Rscript liger.s3s4.R $i
    Rscript liger.s4s3.R $i
    python3 after_liger.py $i
done
