# Combine trait data
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
# cat trait_list.par | awk '{ print $1 }' | while read line; do
    for c in `seq 1 20`; do
        cat data/tmp/$id/$id.$c.dat
    done | awk 'NR == 1 || $1 != "PANEL"' > data/tmp/$id.dat
    cat data/tmp/$id/*.report | awk 'NR == 1 || $1 != "FILE"' > data/tmp/$id.dat.post.report
    cat data/tmp/$id/*.top | awk 'NR == 1 || $1 != "PANEL"' > data/tmp/$id.top.dat
    echo $id
done

# Generate traits info file
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
    avgchisq=`tail -n+2 data/tmp/$id.dat | awk '{ print $19^2 }' | awk -f scripts/avg.awk`
    cat data/tmp/$id/*.report | awk -v chi=$avgchisq -v id=$id 'BEGIN { loc=0; tothit=0; tot=0; } $1 != "FILE" { tothit += $5; tot+=$6; loc++; } END { print id,loc,tot,tothit,chi }'
done | awk 'BEGIN { print "ID NUM.LOCI NUM.JOINT.GENES NUM.GENES AVG.CHISQ" } { print $0 }'  | tr ' ' '\t' > data/traits.par.nfo

# Generate genes info file
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
    tail -n+2 data/tmp/$id.top.dat | cut -f3 | uniq | sort | uniq
done | sort | uniq -c | awk '{ print $2,$1 }' > data/genes.nfo
tail -n+2 data/all.models.par | cut -f3 | sort | uniq -c | awk '{ print $2,$1 }' > data/genes.models.nfo

