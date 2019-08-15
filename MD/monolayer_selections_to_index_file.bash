#!/bin/bash
#make index whose index are not in order

dnmp=$1
resnm=$2
selec_nm=$3
zmiddle=$4
endgrp=$5

#echo 'source monolayer_selections_to_index_file.tcl' > temp.tcl
#echo 'ndxappend md_x1_'${dnmp}'.gro index_'${endgrp}'_'${dnmp}'.ndx '${resnm}' "name '${selec_nm}'" '${zmiddle}'' >> temp.tcl
#echo 'exit' >> temp.tcl

#vmd -dispdev text -e temp.tcl

echo 'source monolayer_selections_to_index_file.tcl' > temp.tcl
echo 'ndxappend md_x1_'${dnmp}'.gro x_nonwater_'${dnmp}'.ndx '${resnm}' "name '${selec_nm}'" '${zmiddle}'' >> temp.tcl
echo 'exit' >> temp.tcl
vmd -dispdev text -e temp.tcl
