#!/bin/bash
 
rel_cutoffs="10 20 30 40 50 60 70 80 90 100" 

geo_file=geometry.xyz
template_file=template.inp
input_file=suppl.inp
 
cutoffs=300
 
for ii in $rel_cutoffs ; do
    work_dir=relCutoff_${ii}Ry
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_rel_cutoff/${ii}/g" \
        -e "s/LT_cutoff/${cutoffs}/g" \
        $template_file > $work_dir/$input_file
    cp $geo_file $work_dir
done
