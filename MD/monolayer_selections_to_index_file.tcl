proc ndxappend {infile ndxfile resnm selection zmiddle} {

    # Load input file (gro, pdb)
    set basemol [mol new $infile]

    # Create atom selections for upper (above zmiddle) and lower monolayers (below zmiddle)
    set sel_upper [atomselect $basemol "resname $resnm and $selection and z > $zmiddle"]
    set sel_lower [atomselect $basemol "resname $resnm and $selection and z < $zmiddle"]
    
    # Numbers of atoms in upper and lower selections
    set nupper [$sel_upper num]
    set nlower [$sel_lower num]

    # Indices for atoms in upper and lower selections
    set index_upper [lsort -integer [$sel_upper list]]
    set index_lower [lsort -integer [$sel_lower list]]
    $sel_upper delete
    $sel_lower delete
    
    # Add 1 to indices since GROMACS indices start at 1
    set i 0
    while {$i < $nupper} {
        set index_upper [lreplace $index_upper $i $i [expr {[lindex $index_upper $i] + 1}]]
        set index_lower [lreplace $index_lower $i $i [expr {[lindex $index_lower $i] + 1}]]
        incr i
    }
    
    # Open index file
    set fid [open $ndxfile a]
    
    # Remove spaces from selection string
    regsub -all {\s} $selection {} selstr

    # Append new index groups to index file for upper and lower selections
    puts $fid "\[ ${resnm}_${selstr}_upper \]"
    puts $fid $index_upper 
    
    puts $fid "\[ ${resnm}_${selstr}_lower \]"
    puts $fid $index_lower
    
    close $fid

}
