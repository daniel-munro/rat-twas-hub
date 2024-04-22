# Combine trait data
# tail -n+2 traits.par | awk '{ print $1 }' | sed 's/.dat//' | while read line; do
cat trait_list.par | awk '{ print $1 }' | while read line; do
    for c in `seq 1 20`; do
        cat tmp/$line/$line.$c.dat
    done | awk 'NR == 1 || $1 != "PANEL"' > tmp/$line.dat
    cat tmp/$line/*.report | awk 'NR == 1 || $1 != "FILE"' > tmp/$line.dat.post.report
    cat tmp/$line/*.top | awk 'NR == 1 || $1 != "PANEL"' > tmp/$line.top.dat
    echo $line
done

# Generate traits info file
tail -n+2 traits.par | while read line; do
    dat=`echo $line | awk '{ print $1 }' | sed 's/.dat//'`
    id=`echo $line | awk '{ print $2 }'`
    avgchisq=`tail -n+2 $dat.dat | awk '{ print $19^2 }' | awk -f avg.awk`
    cat $dat/*.report | awk -v chi=$avgchisq -v id=$id 'BEGIN { loc=0; tothit=0; tot=0; } $1 != "FILE" { tothit += $5; tot+=$6; loc++; } END { print id,loc,tot,tothit,chi }'
done | awk 'BEGIN { print "ID NUM.LOCI NUM.JOINT.GENES NUM.GENES AVG.CHISQ" } { print $0 }'  | tr ' ' '\t' > traits.par.nfo

# Generate genes info file
tail -n+2 traits.par | cut -f1 | sed 's/.dat//' | while read line; do
    tail -n+2 $line.top.dat | cut -f3 | uniq | sort | uniq
done | sort | uniq -c | awk '{ print $2,$1 }' > genes.nfo
tail -n+2 all.models.par | cut -f3 | sort | uniq -c | awk '{ print $2,$1 }' > genes.models.nfo

