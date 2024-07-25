#!/bin/bash

for i in `ls *posterior.txt`; do 
S=`echo $i | cut -f 2 -d "_" | sed 's/s//g'`; 
R=`echo $i | cut -f 3 -d "_" | sed 's/r//g'`;
I=`echo $i | cut -f 4 -d "_"`
V=`awk '$1==101 {print $1,$4}' $i`; 
echo $S $R $I $V; 
done
