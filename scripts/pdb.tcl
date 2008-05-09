#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#  
proc writepsf { file {N_P -1} {MPC -1} {N_CI -1} {N_pS -1} {N_nS -1} } {
    set f [open $file "w"]
    puts $f "PSF"
    puts $f [format "%8d !NATOM" [setmd n_part]]
    # write atoms and create bondlist
    set cnt 1
    set mp [setmd max_part]
    set bondlist {}
    set bondcnt 0
    if { $N_P==-1 } {
	for {set p 0} { $p <= $mp } { incr p } {
	    set tp [part $p p t]
	    if { $tp != "na" } {
		set l [format "%8d T%03d %4d UNX  FE   FE"  $cnt $tp $p]
		if {[string length $l] < 50} {
		    set pad [string repeat " " [expr 50 - [string length $l]]]
		    set l "$l$pad"
		}
		puts $f $l
		set count($p) $cnt
		set bonds [part $p p b]
		foreach bb $bonds {
		    foreach b $bb {  # use 'set b [lindex $bb 0]' instead if you don't wanna have more than the 1st bond
			if {[llength $b ] == 2} { incr bondcnt; lappend bondlist "{$cnt [lindex $b 1]}" }
		    }
		}
		incr cnt
	    }
	}
    } elseif {$N_P=="-molecule"} {
	for {set p 0} { $p <= $mp } { incr p } {
	    set tp [part $p p t]
	    set mid [part $p p mol]
	    if { $tp != "na" } {
		set l [format "%8d T%03d %4d %03d  FE   FE"  $cnt $tp $p $mid]
		if {[string length $l] < 50} {
		    set pad [string repeat " " [expr 50 - [string length $l]]]
		    set l "$l$pad"
		}
		puts $f $l
		set count($p) $cnt
		set bonds [part $p p b]
		foreach bb $bonds {
		    foreach b $bb {  # use 'set b [lindex $bb 0]' instead if you don't wanna have more than the 1st bond
			if {[llength $b ] == 2} { incr bondcnt; lappend bondlist "{$cnt [lindex $b 1]}" }
		    }
		}
	    }
	    incr cnt
	}	
    } else {
	# include topology information in the data stream
	set ids 0
	set j1 [eval list        0         [expr $N_P*$MPC]       [expr $N_P*$MPC+$N_CI]       [expr $N_P*$MPC+$N_CI+$N_pS]       ]
	set j2 [eval list [expr $N_P*$MPC] [expr $N_P*$MPC+$N_CI] [expr $N_P*$MPC+$N_CI+$N_pS] [expr $N_P*$MPC+$N_CI+$N_pS+$N_nS] ]
	foreach ja $j1 je $j2 {
	    for {set p $ja} {$p < $je} {incr p} {
		set tp [part $p p t]
		if { $tp != "na" } {
		    set l [format "%8d T%03d %4d %03d  FE   FE"  $cnt $tp $p $ids]
		    if {[string length $l] < 50} {
			set pad [string repeat " " [expr 50 - [string length $l]]]
			set l "$l$pad"
		    }
		    puts $f $l
		    set count($p) $cnt
		    set bonds [part $p p b]
		    foreach bb $bonds {
			foreach b $bb {  # use 'set b [lindex $bb 0]' instead if you don't wanna have more than the 1st bond
			    if {[llength $b ] == 2} { incr bondcnt; lappend bondlist "{$cnt [lindex $b 1]}" }
                        }
		    }
		}
		incr cnt; if { ($p < $N_P*$MPC-1) && !(($p+1) % $MPC) } { incr ids }
	    }
	    incr ids
	}
    }
    #  write bonds
    puts $f [format "%8d !NBOND                                      " $bondcnt]
    set bondlinecnt 0
    foreach b $bondlist {
	set b [lindex $b 0]
	if {[llength $b] == 2} {
	    incr bondcnt
	    eval "set p1 [lindex $b 0]"
	    eval "set p2 \$count([lindex $b 1])"
	    puts -nonewline $f [format "%8d%8d" $p1 $p2]
	    incr bondlinecnt
	    if {$bondlinecnt == 4} {
		puts $f ""
		set bondlinecnt 0
	    }
	}
    }
    close $f
}

proc writepdb {args} {
    set file [lindex $args 0]
    set de "pos"
    set mode "w"
    set tag ""
    set lscale 1.0
    set namelist ""

    set args [lrange $args 1 end]
    while {$args != ""} {
	set arg [lindex $args 0]
	if {$arg == "-folded"} {
	    set de "folded"
	    set args [lrange $args 1 end]
	} elseif {$arg == "-append"} {
	    set mode "a"
	    set args [lrange $args 1 end]
	} elseif {$arg == "-tag"} {
	    set tag [lindex $args 1]
	    set args [lrange $args 2 end]
	} elseif {$arg == "-lscale"} {
	    set lscale [lindex $args 1]
	    set args [lrange $args 2 end]
	} elseif {$arg == "-names"} {
	    set namelist "[lindex $args 1]"
	    set args [lrange $args 2 end]
	} else {
	    error "unknown flag \"$arg\""
	}
    }
    set f [open $file $mode]
    puts $f "REMARK generated by Espresso"
    if {$tag != ""} {
	puts $f "REMARK $tag"
    }
    # write atoms
    set mp [setmd max_part]
    set cnt 1
    for {set p 0} { $p <= $mp } { incr p } {
	set tp [part $p p t]
	if { $tp != "na" } {
	    set pos [part $p p $de]
	    if { $lscale != 1.0 } {
	       for {set i 0} {$i<3} {incr i} {
	          lset pos $i [expr $lscale*[lindex $pos $i]]
	       }
	    }
	    if { "$namelist" == "" } {
	       set name "FE"
	    } else {
	       set name [lindex $namelist $tp]
	       #names longer than 4 chars will destroy pdb
	       set name [string range $name 0 3]
	    }
	    puts $f [format "ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d" \
			 $cnt $name [expr $p % 10000] [lindex $pos 0] [lindex $pos 1] [lindex $pos 2] $tp]
	    incr cnt
	}
    }
    puts $f "END"
    close $f
}

# This routine writes a pdb file by folding according to molecule information #
proc writepdbfoldtopo { file  {shift 0 } } {

    if { $shift == 0 } {
	set coord [analyze get_folded_positions -molecule ]
    } else {
	set shiftx [lindex $shift 0]
	set shifty [lindex $shift 1]
	set shiftz [lindex $shift 2]
	set coord [analyze get_folded_positions -molecule shift $shiftx $shifty $shiftz ]
    }

    # Get the topology info
    set topo [analyze set]
    set n_molecules [llength $topo]

    # Now write out the particle positions to a pdb file

    set f [open $file "w"]
    puts $f "REMARK generated by Espresso"
    # write atoms
    set mp [setmd max_part]
    set cnt 0

#    set maxparts1 [expr $mp - $nparts2]

    for {set p 0} { $p <= $mp } { incr p } {
	set tp [part $p p t]
	if { $tp != "na" } {
	    set pos [lindex $coord $p]
	    puts $f [format "ATOM %6d  FE  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d" \
			 $cnt [expr $p % 10000] [lindex $pos 1] [lindex $pos 2] [lindex $pos 3] $tp]
	    incr cnt
	}
    }

#    for { set m 0 } { $m < $n_molecules } { incr m } {
#	for { set p 1 } { $p < [llength [lindex $topo $m] ] } { incr p } { 
#	    set pid [lindex [lindex $topo $m] $p]
#	    set tp [part $pid p t]
#	    set pos [lindex $coord [part $pid p id]]
#	    puts $f [format "ATOM %6d  FE  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d" \			 
#    $cnt [expr $cnt % 10000] [lindex $pos 1] [lindex $pos 2] [lindex $pos 3] $tp]
#	    incr cnt
#	}
# }
    close $f
}




proc writepdbfoldchains { file  { chain_start } { n_chains } { chain_length } { box_l } } {
# First collect all the particle coordinates into a tcl list
    set mp [setmd max_part]
    for {set p 0} { $p <= $mp } { incr p } {
	set tp [part $p p t]
	if { $tp != "na" } {
	    lappend coord [part $p p p]
	}
    }

    
    for {set p 0} { $p <= $mp } { incr p } {
	#Identify the chain corresponding to the particle
	set chainid [expr int(floor([expr ($p - $chain_start)/$chain_length]))];

	if { [expr [expr double($p - $chain_start)/double($chain_length)] - $chainid] == 0 } {
	    lappend cm_pos {0 0 0} 
	}
	# Calculate and store the centers of mass for this chain
	for {set j 0} { $j < 3 } { incr j } {
	    lset cm_pos $chainid $j [expr [ lindex $cm_pos $chainid $j ] + double([lindex $coord $p $j])/double($chain_length)]
	}
    }

# Now fold the particles by chain
    for { set chainid 0 } { $chainid < $n_chains } { incr chainid } {
	for { set j 0 } { $j < 3 } { incr j } {
	    set Lj [lindex $box_l $j]
	    set arg [expr double([lindex $cm_pos $chainid $j])/$Lj ]
	    set tmp [expr floor($arg)]
	    set cm_tmp 0.0
	    for { set i 0 } { $i < $chain_length } { incr i } {
		set index [expr $chainid*$chain_length + $i] 
		lset coord $index $j [expr [lindex $coord $index $j] - $tmp*$Lj]
		set $cm_tmp  [expr $cm_tmp + [lindex $coord $index $j]/$chain_length];
	    }
	    if { $cm_tmp < 0 || $cm_tmp > $Lj } {
		put "chain center of mass is out of range"
	    }
	}
    }

# Now write out the particle positions to a pdb file

    set f [open $file "w"]
    puts $f "REMARK generated by Espresso"
    # write atoms
    set mp [setmd max_part]
    set cnt 0
    for {set p 0} { $p <= $mp } { incr p } {
	set tp [part $p p t]
	if { $tp != "na" } {
	    set pos [lindex $coord $p]
	    puts $f [format "ATOM %6d  FE  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d" \
			 $cnt [expr $p % 10000] [lindex $pos 0] [lindex $pos 1] [lindex $pos 2] $tp]
	    incr cnt
	}
    }
    close $f
}

