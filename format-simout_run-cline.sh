#!/bin/bash

# Make output directories
mkdir vcfs
mkdir genos
mkdir formatted_genos
mkdir posteriors

# Format the genos files
for i in `ls *vcf`; do
NAME=`echo $i | sed 's/\.vcf//g'`
grep -v "##" $i | cut -f1,2,10- | sed 's/#//g' > ${NAME}'_genos.txt'
done

# Run gghybrid cline analysis
for i in `ls *genos.txt`; do
NAME=`echo $i | sed 's/_genos.txt//g'`
Rscript ~/other_projs/chickadee_hybridzone_movement/gghybrid_chickadee_moving.R $NAME >> $1 2>&1
mv $NAME'.vcf' vcfs/
mv $NAME'_genos.txt' genos/
mv $NAME'_genos_formatted.txt' formatted_genos/
mv $NAME'_clines_posterior.txt' posteriors/
done
