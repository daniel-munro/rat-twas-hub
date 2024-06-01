# Merge trait data and generate summary files

# Inputs:
# data/traits.par
# data/twas_out/{trait}/{trait}.{1..20}.dat
# data/twas_out/{trait}/{trait}.{1..20}.top
# data/twas_out/{trait}/{trait}.{1..20}.post.report

# Outputs:
# data/twas_out/{trait}.dat
# data/twas_out/{trait}.dat.post.report
# data/twas_out/{trait}.top.dat
# data/traits.par.nfo
# data/genes_n_hits.nfo
# data/genes_n_models.nfo

# Combine trait data
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
    for c in `seq 1 20`; do
        cat data/twas_out/$id/$id.$c.dat
    done | awk 'NR == 1 || $1 != "PANEL"' > data/twas_out/$id.dat
    cat data/twas_out/$id/*.report | awk 'NR == 1 || $1 != "FILE"' > data/twas_out/$id.dat.post.report
    cat data/twas_out/$id/*.top | awk 'NR == 1 || $1 != "PANEL"' > data/twas_out/$id.top.dat
    echo $id
done

# Generate traits info file
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
    avgchisq=`tail -n+2 data/twas_out/$id.dat | awk '{ print $19^2 }' | awk 'BEGIN { n=0; s=0; } { n++; s += $1; } END { if ( n > 0 ) print s / n; else print "NA"; }'`
    cat data/twas_out/$id/*.report | awk -v chi=$avgchisq -v id=$id 'BEGIN { loc=0; tothit=0; tot=0; } $1 != "FILE" { tothit += $5; tot+=$6; loc++; } END { print id,loc,tot,tothit,chi }'
done | awk 'BEGIN { print "ID NUM.LOCI NUM.JOINT.GENES NUM.GENES AVG.CHISQ" } { print $0 }' | tr ' ' '\t' > data/traits.par.nfo

# Generate genes info files
tail -n+2 data/traits.par | awk '{ print $2 }' | while read id; do
    tail -n+2 data/twas_out/$id.top.dat | cut -f3 | grep -o '^[^:.]*' | uniq | sort | uniq
done | sort | uniq -c | awk '{ print $2,$1 }' > data/genes_n_assoc.nfo
tail -n+2 data/all_models.par | cut -f3 | grep -o '^[^:.]*' | sort | uniq -c | awk '{ print $2,$1 }' > data/genes_n_models.nfo

