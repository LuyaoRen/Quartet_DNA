cat benchmark.men.vote.diffbed.filtered | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$7"\t.\t.\t.\tGT\t"$6}' | grep -v '2_y' > LCL5.body
cat benchmark.men.vote.diffbed.filtered | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$15"\t.\t.\t.\tGT\t"$14}' | grep -v '2_y' > LCL6.body
cat benchmark.men.vote.diffbed.filtered | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$23"\t.\t.\t.\tGT\t"$22}' | grep -v '2_y' > LCL7.body
cat benchmark.men.vote.diffbed.filtered | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$31"\t.\t.\t.\tGT\t"$30}' | grep -v '2_y' > LCL8.body

for i in *txt; do cat $i | awk '{ if ((length($3) == 1) && (length($4) == 1)) { print } }' | grep -v '#' | cut -f3,4 | sort |uniq -c | sed 's/\s\+/\t/g' | cut -f2 > $i.mut; done