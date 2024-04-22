cat trait_list.par | while read line; do
    g=`echo $line | awk '{ print $1 }'`
    n=`echo $line | awk '{ print $2 }'`
    sbatch -n 1 -t 04:00:00 --mem-per-cpu=8G --array=1-20 --job-name="h2" TWAS.sh $g $n
done
