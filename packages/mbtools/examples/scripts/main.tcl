#-----------------------------------------------------------#
#
# Simple example script for running a membrane simulation
#
# Implements settings contained in a parameter file which should be
# given as an argument.
#
# -----------------------------------------------------------#

# Get the name of the current namespace
set this [namespace current]


# --- ensure required packages are available  ---------------#
set result [package require ::mmsg]
 ::mmsg::send $this "loaded version [package require ::mbtools::analysis] of analysis"
 ::mmsg::send $this "loaded version [package require cmdline] of cmdline"
 ::mmsg::send $this "loaded version [package require ::mbtools::system_generation] of system_generation"

# Capture the child namespaces of analysis and system_generation so
# that we can explicitly allow messages from these.
set message_allowlist { :: ::mbtools::utils ::mbtools::system_generation ::mbtools::analysis }
set children [namespace children ::mbtools::analysis]
foreach child $children {
    lappend message_allowlist $child
}
set children [namespace children ::mbtools::system_generation]
foreach child $children {
    lappend message_allowlist $child
}

# Set the namespaces from which messages will be printed
catch { ::mmsg::setnamespaces $message_allowlist }




# ---- Process Command Line Args --------------------------- #

# Note that this can cause problems reading checkpoints.  The usage
# below seems to work ok.
set options {
    {n.arg      1    set the number of processors }
}
set usage "Usage: main.tcl n: paramsfile"
array set params [::cmdline::getoptions argv $options $usage]
set paramsfile [lindex $argv 0]

# Enable debugging messages
#::mmsg::enable debug

#----------- System Parameters ----------------------------#
# Set a bunch of default parameters.
set thermo Langevin
set warmup_temp 0
set warmsteps 100
set warmtimes 20
set free_warmsteps 0 
set free_warmtimes 0 
set startj 0    
set startk 0
set startmdtime 0
set npt off
set mgrid 8
set stray_cut_off 1000.0
set use_vmd "offline"
set tablenames ""
set trappedmols ""
::mmsg::send $this "using paramsfile: $paramsfile"
source $paramsfile



#----------- Default Parameters set from System Params ----#
if { $warmup_temp == 0 } {
    set warmup_temp [lindex $systemtemp 0 ]
}



# ----------- Initialization ------------------ -----------#

# Attempt to read a checkpoint file
set checkpointexists [ ::mbtools::utils::readcheckpoint $outputdir ]

# Set the starting time for this run ie override the value in checkpoint file
set starttime [clock seconds]

if { !$checkpointexists } {
    # No checkpoint exists so we need to setup everything from scratch

    # Setup the output directory by creating it and copying
    # forcetables to it
    ::mbtools::utils::setup_outputdir  $outputdir -paramsfile $paramsfile -tabdir $tabledir -tabnames $tablenames
	
    # Set the box dimensions
    setmd box_l [lindex $setbox_l 0] [lindex $setbox_l 1] [lindex $setbox_l 2]
	
    # Specify the bonded interactions
    ::mbtools::utils::set_bonded_interactions $bonded_parms

    # Specifiy the non-bonded tabulated interactions .. actually this
    # routine is only kept for backwards compatibility other_nb should
    # be used instead

#    ::setup_utilities::set_tab_interactions $tablenames $tabledir $outputdir 

    # Specify any other non-bonded interactions
    if { [ catch { ::mbtools::utils::set_nb_interactions $nb_interactions } ] } {
	mmsg::send $this "no non-bonded interactions used"
    }


    set topology [::mbtools::system_generation::setup_system $system_specs $setbox_l $moltypes]
    # Set the generated topology into the internals of espresso.
    ::mbtools::utils::set_topology $topology
    set trappedmols [::mbtools::system_generation::get_trappedmols]
    # Fix molecules if necessary
    if { $trappedmols != -1 } {
	::mbtools::utils::trap_mols $trappedmols
    }

    #Initialise Random Number Generator
    #----------------------------------------------------------#
    ::mbtools::utils::init_random $params(n)

    # ----------- Integration Parameters before warmup -----------#
    setmd time_step $warm_time_step
    setmd skin      $verlet_skin
    setmd gamma     $langevin_gamma
    setmd temp      $warmup_temp
    
    # Set the topology and molecule information
    #----------------------------------------------------------#
    

    #write topology file
    set f_topo [open "$outputdir/$ident.top" w]
    blockfile_write_topology $f_topo write topology   
    close $f_topo

   
    # Check if there are any extra vmdcommands and if not initialize a default
    if { [catch { set dummy $vmdcommands } ] } {
	set vmdcommands ""
    }
    ::mbtools::utils::initialize_vmd $use_vmd $outputdir $ident -extracommands $vmdcommands

    #Perform the warm up integration
    #----------------------------------------------------------#
    mmsg::send $this "warming up at [setmd temp]"

    

    ::mbtools::utils::warmup  $warmsteps $warmtimes -cfgs 10 -outputdir $outputdir

    # Bilayer systems may have been setup with fixed z positions for
    # particles.  Here we make sure that all particles are unfixed
    # after warmup
    set userfixedparts [::mbtools::system_generation::get_userfixedparts ]
    for {set i 0} { $i <  [setmd n_part] } {incr i} {
	if { [lsearch $userfixedparts $i ] == -1 } {
	    part [expr $i] fix 0 0 0
	}
    }

    setmd time_step $main_time_step
    setmd temp     [lindex $systemtemp 0]
    ::mmsg::send $this "warming up again at  [setmd temp]"
    ::mbtools::utils::warmup $free_warmsteps $free_warmtimes -startcap 1000 -outputdir $outputdir
    

    # ----------- Integration Parameters after warmup -----------#
    setmd time_step $main_time_step
    setmd temp     [lindex $systemtemp 0]
 
    ::mbtools::analysis::setup_analysis $analysis_flags -outputdir  $outputdir -g $mgrid -str $stray_cut_off
    
    mmsg::send $this "starting integration: run $int_n_times times $int_steps steps"

    # Reset the time to a starttime (usually zero) after warmup
    setmd time $startmdtime   

} else {

    # A checkpoint exists so all we need to do is reset the topology and setup analysis again
    ::mbtools::utils::read_topology "$outputdir/$topofile"
    ::mbtools::analysis::setup_analysis $analysis_flags -outputdir  $outputdir -g $mgrid -str $stray_cut_off
    ::mbtools::utils::initialize_vmd $use_vmd $outputdir $ident

    # Yikes I hope this is right.  We want to make sure that we start
    # exactly from where the checkpoint was written
    set startj $j
    set startk [expr $k + 1]

}

#Main Integration                                          #
#----------------------------------------------------------#

set j $startj

setmd temp $systemtemp
if { $thermo == "DPD" } {
    thermostat off
    thermostat set dpd $systemtemp $dpd_gamma $dpd_r_cut
    mmsg::send $this "DPD thermostat has been set"
    mmsg::send $this "Thermostat is: [thermostat]"
}
if { $npt == "on" } {
    integrate set npt_isotropic $p_ext $piston_mass 1 1 0
    mmsg::send $this "npt integrator has been set"
    flush stdout
    #-cubic_box
    thermostat set npt_isotropic $systemtemp  $gamma_0  $gamma_v
}

set timingstart [clock clicks -milliseconds]
for {set k $startk } { $k <  $int_n_times } { incr k} {

    mmsg::send $this "run $k at time=[setmd time]"



    # Call all of the analyze routines that we specified when setting
    # up our analysis
    ::mbtools::analysis::do_analysis

    # If k is a multiple of analysis_write_frequency then write the
    # analysis results to file
    if { $k%$analysis_write_frequency ==0 } {
	::mbtools::analysis::print_averages
    }

    # If k is a multiple of write_frequency then write out a full
    # particle configuration
    if { $k%$write_frequency==0 } {
	polyBlockWrite "$outputdir/$ident.[format %04d $j].gz" {time box_l npt_p_diff } {id pos type v f molecule} 
	mmsg::send $this "wrote file $outputdir/$ident.[format %04d $j].gz " 
	flush stdout

	if { $use_vmd == "offline" } {
	    writepdbfoldtopo "$outputdir/$ident.vmd[format %04d $j].pdb"  
	}

	incr j

    }


    # Do the real work of integrating equations of motion
    integrate $int_steps
    set timingcurr [clock clicks -milliseconds]

    set elapsedtime [expr  $timingcurr - $timingstart]
    ::mmsg::send $this "elapsed time: $elapsedtime"

    # Write a checkpoint to allow restarting.  Overwrites previous
    # checkpoint
    mmsg::send $this "setting checkpoint $k [setmd time] $j"    
    checkpoint_set "$outputdir/checkpoint.latest.gz"

    # Try to copy a checkpoint to the backup checkpoint folder.
    # Usefull if the program crashes while writing a checkpoint
    if { [ catch { exec cp "$outputdir/checkpoint.latest.gz" "$outputdir/checkpoint_bak/checkpoint.latest.gz" } ] } {
	mmsg::warn $this "warning: couldn't copy backup checkpoint"
    }

}

# terminate program
mmsg::send $this "\n\nfinished"
exit 1



