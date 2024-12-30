### Calculate the h2 of target set
#!/bin/bash
for file in *snplist; do
    cd "/home/zdzhao/mulit_omics_result/${file%.snplist}"
    echo "${file%.snplist}" "no_${file%.snplist}" > "${file%.snplist}".mgrm
    gcta64 --bfile /home/zdzhao/1577_cattle/1577 --extract "$file" --make-grm --out "${file%.snplist}"
    gcta64 --bfile /home/zdzhao/1577_cattle/1577 --exclude "$file" --make-grm --out "no_${file%.snplist}"
    gcta64 --mgrm "${file%.snplist}".mgrm --reml --pheno /home/zdzhao/phenotype/ph.txt --mpheno $i --out "tr$i" --reml-alg 2 --threads 40 --reml-maxit 10000
done

### Demo of preparing a document
for file in tr*.hsq; do
    sed -n '8,11p' "$file" | awk -v filename="$file" '{print filename, \$0 >> "ld.hsq"}'; done
sed -i 's/\s\+/\t/g' ld.hsq 

##### Change the names
awk 'BEGIN { FS=OFS="\t" } { 
    if (\$2 == "V(G1)/Vp") \$2 = "LD1"
    else if (\$2 == "V(G2)/Vp") \$2 = "LD2"
    else if (\$2 == "V(G3)/Vp") \$2 = "LD3"
    else if (\$2 == "V(G4)/Vp") \$2 = "LD4"
    print
}' ld.hsq >ld.h2_sum

### Calculate the number of target set
for file in *snplist; do
    size=$(wc -l < "$file")
    echo -e "${file%.snplist}\t$size"
done > demo.size.grm
### Calculate the variants score
Rscript score.R

### Selection of variants based on thresholds
for i in $(seq 5 10 30); do
    sort -k2 -n demo_score.txt > sorted_scores.txt
    line=$(wc -l < sorted_scores.txt)
    line_bottom_i=$((i * line / 100))
    line_top_i=$(((100 - i) * line / 100))
    awk -v line_i="$line_bottom_i" 'NR < line_i {print \$1}' sorted_scores.txt > score_bottom_$i.snplist
    awk -v line_i="$line_top_i" -v line="$line" 'NR > line_i {print \$1}' sorted_scores.txt > score_top_$i.snplist
done
