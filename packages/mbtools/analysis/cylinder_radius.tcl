# ::mbtools::analysis::analyze_cylinder_radius --
#
# Calculate the radius of a cylinder
#   

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::cylinder_radius {
    variable R_middle 0
    variable R_inner 0
    variable R_outer 0

    variable verbose

    namespace export printav_cylinder_radius
    namespace export setup_cylinder_radius
    namespace export analyze_cylinder_radius
    namespace export resetav_cylinder_radius
}

proc ::mbtools::analysis::cylinder_radius::resetav_cylinder_radius { } {
    variable R_middle 0
    variable R_inner 0
    variable R_outer 0
}

proc ::mbtools::analysis::cylinder_radius::printav_cylinder_radius { } {
    variable R_middle
    variable R_inner
    variable R_outer
    variable f_tvsen
    global ::mbtools::analysis::time
    puts -nonewline $f_tvsen "$time $R_middle $R_inner $R_outer"
    puts $f_tvsen ""
    flush $f_tvsen
}

proc ::mbtools::analysis::cylinder_radius::setup_cylinder_radius { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsen
    variable verbose

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_cylinder_radius verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_cylinder_radius$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    ::mmsg::debug  [namespace current]  "opening $outputdir/time_vs_cylinder_radius$suffix "
    set f_tvsen [open "$outputdir/time_vs_cylinder_radius$suffix" $iotype ]

    if { $newfile || $iotype == "w"} {
	puts $f_tvsen "\# Components of the cylinder_radius "
	puts $f_tvsen "\# Time R_middle R_inner R_outer"
    }
    flush $f_tvsen
}


proc ::mbtools::analysis::cylinder_radius::analyze_cylinder_radius {  } {
    ::mmsg::send [namespace current] "analyzing cylinder_radius"
    variable R_middle
    variable R_inner
    variable R_outer
    variable verbose

    set box [setmd box_l]

    set topology [analyze set]

    set com {0. 0.}
    set num_mol 0
    # Assume cylinder is built along the z-direction
    # Calculate center of mass of the cylinder in the x-y plane.    
    # Analysis does not discriminate stray beads.
    foreach mol $topology {
	set headp [lindex $mol 1]
	set tb [lindex $mol end]
	
	set hpos [foldpart [part $headp print pos] $box]
	
	lset com 0 [expr [lindex $com 0]+[lindex $hpos 0]]
	lset com 1 [expr [lindex $com 1]+[lindex $hpos 1]]
	incr num_mol
    }
    lset com 0 [expr [lindex $com 0]/(1.*$num_mol)]
    lset com 1 [expr [lindex $com 1]/(1.*$num_mol)]

    set num_outer 0
    # Now calculate two radii: inner and outer by looking at the orientation of a lipid
    foreach mol $topology {
	set headp [lindex $mol 1]
	set tb [lindex $mol end]

	set hpos [foldpart [part $headp print pos] $box]
	set tpos [foldpart [part $tb    print pos] $box]

	# Define vector CoM to head bead
	set com_h_vec {0. 0.}
	lset com_h_vec 0 [expr [lindex $hpos 0]-[lindex $com 0]]
	lset com_h_vec 1 [expr [lindex $hpos 1]-[lindex $com 1]]

	# Define vector tail bead to head bead
	set t_h_vec {0. 0.}
	lset t_h_vec 0 [expr [lindex $hpos 0]-[lindex $tpos 0]]
	lset t_h_vec 1 [expr [lindex $hpos 1]-[lindex $tpos 1]]

	# dot product between the two vectors
	set com_h_dot_t_h [expr [lindex $com_h_vec 0]*[lindex $t_h_vec 0] + \
			       [lindex $com_h_vec 1]*[lindex $t_h_vec 1]]
	
	if {$com_h_dot_t_h > 0} {
	    incr num_outer
	    # same orientation -- it's the outer leaflet
	    set R_outer [expr $R_outer + sqrt(pow([lindex $hpos 0]-[lindex $com 0],2) + \
						  pow([lindex $hpos 1]-[lindex $com 1],2))]
	} else {
	    # antiparallel orientation -- it's the inner leaflet
	    set R_inner [expr $R_inner + sqrt(pow([lindex $hpos 0]-[lindex $com 0],2) + \
						  pow([lindex $hpos 1]-[lindex $com 1],2))]
	}
    }
	
    set R_outer [expr $R_outer/(1.*$num_outer)]
    set R_inner [expr $R_inner/(1.*($num_mol-$num_outer))]
    set R_middle [expr ($R_inner+$R_outer)/2.]

    if { $verbose } {
	::mmsg::send [namespace current] "cylinder_radius: R_middle: $R_middle -- R_inner: $R_inner -- R_outer: $R_outer"
    }
    
    ::mmsg::debug [namespace current] "done"

}


