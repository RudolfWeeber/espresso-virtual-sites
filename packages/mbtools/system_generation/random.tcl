# 
# Routines for creating a random "gas" of lipids
#

namespace eval ::mbtools::system_generation {}

namespace eval ::mbtools::system_generation::random {
    namespace export create_random
}

# ::system::generation::random::create_random --
#
# puts lipids in random positions with random orientation
#
proc ::mbtools::system_generation::random::create_random { args } {
    ::mmsg::send [namespace current] "placing lipids in a random fluid "
    set options {
	{bondl.arg     1.0   "bond length between atoms"  }
	{shuffle "shuffle topology before placement"}
	{exclude.arg "" "a region where no lipids should be placed"}
    }
    set usage "Usage: create_random \[bondl:shuffle:exclude]"
    array set params [::cmdline::getoptions args $options $usage]
    
    
    global ::mbtools::system_generation::topology
    global ::mbtools::system_generation::boxl
    
    # First work out how many mol types there are and construct a list
    # of their lengths
    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]

    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    set maxtries 1000

    if { $params(shuffle) } {
        set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }
    
    foreach mol $topology {



	set tries 0
	set tailpos { 0 0 0 }
	set isallowed 0
	while { !$isallowed  } {

	    if {  ($tries > $maxtries) } {
		mmsg::err [namespace current] "could not place molecule exceeded max number of tries"
	    }


	    # First we choose a random point in space for the tail atom
	    lset tailpos 0 [expr $bx*[t_random]]
	    lset tailpos 1 [expr $by*[t_random]]
	    lset tailpos 2 [expr $bz*[t_random]]
	    if { [llength $params(exclude)] != 0 } {		
		set isallowed [::mbtools::utils::isoutside $tailpos $params(exclude) ]
	    } else {
		set isallowed 1
	    }
	    incr tries

	}

	# Now choose a random orientation vector.  
	lappend orient [expr 2*[t_random]-1]
	lappend orient [expr 2*[t_random]-1]
	lappend orient [expr 2*[t_random]-1] 

	
	::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
	
	unset tailpos
	unset orient
    }


    # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }

    return

}


























    
