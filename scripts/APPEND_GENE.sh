DAT=$1
NAME=$2
NORM=$3

cat $DAT \
| tail -n+2 \
| awk '{ print NR,$3,$19 }' \
| sort -k2,2 -k1,1g \
| awk -v n=$NORM -v link="$NAME" 'BEGIN {prev=""}
{
    if($2 != prev) {
        if(prev!="") {
            cur_str = link" | "sprintf("%2.1f",sum/tot/n)" | "sprintf("%2.1f",sum/tot)" | "sprintf("%2.1f",max)" "cur_str
            print cur_str >> "jekyll/genes/"prev".md"
        }
        max=0
        sum=0
        tot=0
        cur_str = ""
    }
    cur_str = cur_str" | "sprintf("%2.1f",$3)
    if(max < $3^2) max=$3^2
    sum += $3^2
    tot++
    prev=$2
}
END {
    cur_str = link" | "sprintf("%2.1f",sum/tot/n)" | "sprintf("%2.1f",sum/tot)" | "sprintf("%2.1f",max)" "cur_str
    print cur_str >> "jekyll/genes/"prev".md"
}'
