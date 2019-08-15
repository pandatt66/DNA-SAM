#!/bin/bash

dnmp=$1
input_str='del 1-4\n'
input_str+='del 6-8\n'
input_str+='name 5 dNMP\n'
input_str+='a SH\n'
if [ "$dnmp" = "dCMP" ]; then 
input_str+='a C4 and a C5 and a C6 and a N1 and a C2 and a N3\n'
input_str+='a N1\n'
input_str+='a C4\n'
elif [ "$dnmp" = "dTMP" ]; then
input_str+='a N3 and a C4 and a C5 and a C6 and a N1 and a C2\n'
input_str+='a N1\n'
input_str+='a C4\n'
elif [ "$dnmp" = "dAMP" ]; then
input_str+='a C5 and a C6 and a N1 and a C2 and a N3 and a C4 and a C8 and a N7 and a N9\n'
input_str+='a N1\n'
input_str+='a C8\n'
else
input_str+='a C8 and a N9 and a C4 and a C5 and a N7 and a C6 and a N1 and a C2 and a N3\n'
input_str+='a N1\n'
input_str+='a C8\n'
fi
input_str+='name 7 dNMP_base\n'
input_str+='a P\n'
input_str+='a C13\n'

input_str+='q\n'

gmx make_ndx  -f x_nonwater_${dnmp}.tpr -o x_nonwater_${dnmp}.ndx << EOF

`echo -e ${input_str}`
EOF
