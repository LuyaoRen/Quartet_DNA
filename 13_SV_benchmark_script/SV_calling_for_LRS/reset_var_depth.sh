#!/bin/bash


if [ $1 = "sniffles" ] ; then

#sniffles/cutesv
#perl -lane 'my $depth=10;if(/^#/){print $_}{@format=split /:/,$F[9]}if(($format[1]+$format[2])>=$depth){print $_}' $2>$3

perl -lane 'my $RE=3;if(/^#/){print $_}elsif(/RE=(\d+)/){if($1>=$RE){print $_}}else{next}' $2>$3

elif [ $1 = "svim" ] ; then 

#svim
#perl -lane 'my $depth=10;if(/^#/){print $_}{@format=split /:/,$F[9]}{@allele=split /,/,$format[2]}if($format[1]>=$depth){print $_}' $2>$3

perl -lane 'my $SUPPORT=3;if(/^#/){print $_}elsif(/SUPPORT=([-|\d]+)/){if($1>=$SUPPORT){print}}else{next}' $2>$3

elif [ $1 = "pbsv" ] ; then 
#pbsv
#perl -lane 'my $depth=5;if(/^#/){print $_}{@format=split /:/,$F[9]}{@allele=split /,/,$format[1]}if($allele[1]>=$depth){print $_}' $2>$3

perl -lane 'my $var=3;if(/^#/){print $_}else{@format=split /:/,$F[9];@allele=split /,/,$format[1];if($allele[1]>=$var){print $_}}' $2>$3

elif [ $1 = "nanosv" ] ; then 
#nanosv
#perl -lane 'my $depth=5;if(/^#/){print $_}{@format=split /:/,$F[9]}{@allele=split /,/,$format[1]}if($allele[1]>=$depth){print $_}' $2>$3

perl -lane 'my $var=3;if(/^#/){print $_}else{@format=split /:/,$F[9];@allele=split /,/,$format[2];if($allele[1]>=$var){print $_}}' $2>$3

else

break

fi
