#!/bin/bash
dnmp=$1
time_file=$2
trs_dist=$3
endgrp=$4

#this bash file aims at i)dump all the selectd time for configuration; ii)translate system to put solvent in center of box ;iii)generate tpr files
tms=(`cat ${time_file}|awk '{print $1}' `)
dis=(`cat ${time_file}|awk '{printf ("%.5f\n",$2)}'`)
ctline=`cat ${time_file}|awk '{print NR}' |tail -1`
index=`expr $ctline - 1`

for i in `seq 0 $index`; do
  echo 0| gmx trjconv -f ${dnmp}_ctr_SH.trr -s md_x1_${dnmp}.tpr -dump ${tms[i]} -o pmf_${dnmp}_${dis[i]}.gro
done

for i in `seq 0 $index`; do
  echo 0|gmx trjconv -f pmf_${dnmp}_${dis[i]}.gro -s md_x1_${dnmp}.tpr -pbc mol -n index_${endgrp}_${dnmp}.ndx -trans 0 0 ${trs_dist} -o ${dnmp}_ctr_SOL_${dis[i]}.gro
done

if [ "$endgrp" = "cooc" ]; then 
echo 0| gmx trjconv -f ref_md_boxz_${dnmp}.pdb -s md_x1_${dnmp}.tpr -pbc mol -trans 0 0 ${trs_dist} -o ref_pmf_${dnmp}.pdb
else
echo 0| gmx trjconv -f ref_throw_${dnmp}.pdb -s md_x1_${dnmp}.tpr -pbc mol -trans 0 0 ${trs_dist} -o ref_pmf_${dnmp}.pdb
fi

for i in `seq 0 $index`; do
    grompp -f pull.mdp -c ${dnmp}_ctr_SOL_${dis[i]}.gro -p ${endgrp}_md_${dnmp}.top -r ref_pmf_${dnmp}.pdb -n index_${endgrp}_${dnmp}.ndx -o ${dnmp}_x1_${dis[i]}.tpr -maxwarn 1
done
