# N_TRAITS=`wc -l traits.par | awk '{ print $1 - 1 }'`
# sbatch -n 1 -t 08:00:00 --mem-per-cpu=8G --job-name="report" --array=1-$N_TRAITS REPORT.sh
N_TRAITS=`wc -l traits.par | awk '{ print $1 - 1 }'`
for i in `seq 1 $N_TRAITS`; do
    Rscript REPORT_SINGLE.R $i
done
