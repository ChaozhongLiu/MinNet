for i in $(seq 1 8)
do
    python3 glue.s3s4.py $i
    python3 glue.s4s3.py $i
done
