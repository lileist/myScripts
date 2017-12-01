#!/bin/sh
#This code is used to prepare 'jobfile.sh' script used to run exafs-assisted basin hopping jobs via TACC Launcher
#start=$(printf -v int '%d' "$1")
#end=$(printf -v int '%d' "$2")
#for i in {$start..$end}
for i in {6..11}
do
    run=$(printf "run-%i" $i)
    mkdir=$(printf "mkdir -p run-%i" $i)
    mkdir_local=$(printf "mkdir -p /state/partition1/run-%i" $i)
    cp_pos=$(printf "cp ../random_stru/%i.xyz %s/pos.xyz" "$i" "$run")
    cp_files=$(printf "cp bh.py basin.py PdAu_sc.eam.fs exp_exafs.chi %s" $run)
    run_job=$(printf "cd %s; python bh.py >/state/partition1/$run/out; cd .." $run)
    cmd=$(printf "%s; %s; %s; %s; %s" "$mkdir" "$cp_pos" "$cp_files" "$mkdir_local" "$run_job")
    #cmd=$(printf "mkdir -p run-%i; cp ../random_stru/%i.xyz run-%i/pos.xyz; cp bh.py basin.py PdAu_sc.eam.fs exp_exafs.chi run-%i; cd run-%i; python bh.py >out; cd .." $i)
    echo $cmd >>jobfile_6.sh
done
