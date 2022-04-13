for i in $(seq 1 8)
do
    Rscript liger_online.R $i
    python3 after_ligerOnline.py $i
done
