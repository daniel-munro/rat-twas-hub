# For one trait, for each gene, write summary stats and Z-score per model to the table in that gene's markdown file.
# Output as HTML table rows so that they can be inserted into the table in the HTML layout for the gene page.

# Inputs:
# data/twas_out/{trait}.all.tsv

# Outputs:
# jekyll/genes/{gene}.md (append)

OUT=$1
TRAIT_LINK=$2
AVG_CHISQ=$3

cat $OUT \
| tail -n+2 \
| awk '{ print NR,$3,$19 }' \
| awk '{ $2 = gensub(/[:.].*$/, "", "g", $2); print }' \
| sort -k2,2 -k1,1g \
| awk -v avg_chisq=$AVG_CHISQ -v trait_link="$TRAIT_LINK" 'BEGIN {prev=""}
{
    if($2 != prev) {
        if(prev!="") {
            cur_str = "<tr><td>"trait_link"</td><td>"sprintf("%2.1f",sum/tot/avg_chisq)"</td><td>"sprintf("%2.1f",sum/tot)"</td><td>"sprintf("%2.1f",max)"</td>"cur_str"</tr>"
            print cur_str >> "jekyll/genes/"prev".md"
        }
        max=0
        sum=0
        tot=0
        cur_str = ""
    }
    cur_str = cur_str"<td>"sprintf("%2.1f",$3)"</td>"
    if(max < $3^2) max=$3^2
    sum += $3^2
    tot++
    prev=$2
}
END {
    cur_str = "<tr><td>"trait_link"</td><td>"sprintf("%2.1f",sum/tot/avg_chisq)"</td><td>"sprintf("%2.1f",sum/tot)"</td><td>"sprintf("%2.1f",max)"</td>"cur_str"</tr>"
    print cur_str >> "jekyll/genes/"prev".md"
}'
