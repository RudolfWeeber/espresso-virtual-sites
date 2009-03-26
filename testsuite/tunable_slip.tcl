# This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
# It is therefore subject to the ESPResSo license agreement which you
# accepted upon receiving the distribution and by which you are
# legally bound while utilizing this file in any form or way. 
# There is NO WARRANTY, not even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# You should have received a copy of that license along with this
# program; if not, refer to http://www.espresso.mpg.de/license.html
# where its current version can be found, or write to
# Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 
# 55021 Mainz, Germany.
# Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#############################################################
#                                                           #
# Basic tests of the tunable-slip boundary interactions     #
#                                                           #
# 1) check constraints                                      #
# 2) check implementation with DPD                          #
#                                                           #
# Reference:                                                #
# J.Smiatek, M.P. Allen, F. Schmid:                         #
# "Tunable-slip boundaries for coarse-grained simulations   #
# of fluid flow", Europ. Phys. J. E 26, 115 (2008)          #
#                                                           #
# Responsible for the testcase and the source code:         #
# Jens Smiatek, smiatek@physik.uni-bielefeld.de             #
#                                                           #
#############################################################

set errf [lindex $argv 1]

proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
        set f [open $errf "w"]
        puts $f "not compiled in: $feature"
        close $f
        exit -42
    }
}

puts "----------------------------------------"
puts "- Testcase tunable_slip.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

set errf [lindex $argv 1]

require_feature "TUNABLE_SLIP"
require feature "DPD"

# System parameters 
set box_l                10.0

set box_x                $box_l

set box_y                $box_l

set box_z                $box_l

# reuse of Verlet lists  
setmd verlet_reuse

# skin depth
setmd skin               0.4

# box length
setmd box_l              $box_l $box_l $box_l

# volume of system
set volume               [expr $box_x*$box_y*$box_z]

# periodic boundary conditions
setmd periodic           1 1 1

# cell size
setmd cell_size

# cell grid
setmd cell_grid  

################################################
#         Solvent parameters                   #
################################################

# type number for solvent molecules
set solvent_id           0

# density for solvent molecules
set density_s            3.0

# number of solvent particles
set n_solvent            [expr int (ceil($density_s*$volume))]

# total number of particles
set n_part               $n_solvent

################################################
#         Interaction parameters               #
################################################

# Lennard-Jones ################################

# epsilom
set lj_eps               1.0

# sigma
set lj_sig               1.0

# cutoff
set lj_cut               1.12246

# shift
set lj_shift             0.0

# offset
set lj_off               0.0

################################################
#            Thermostat                        #
################################################
# temperature
setmd temperature

set temp                 1.0

#Tunable Slip Boundaries
set gamma_L              1.0

# cutoff has to be larger than LJ-interaction range
set r_cut_L              2.0

# Extrem accelaration of the flow
set vx                   10.0

set vy                   0.0

set vz                   0.0

# DPD thermostat ###############################

set dpd_temperature            1.0

# DPD friction coefficient
set dpd_gamma                  5.0

# DPD-r_cutoff radius
set dpd_r_cut                  1.0

# length of time step
setmd time_step                 0.01

set timestep                    0.01

# Length of integration scheme
set int_steps                   100

# Length of integration run
set int_loops                   5

# precision of data output
set tcl_precision     6

# initialize random number generator
set t_random_seed     [expr int(rand()*99999)^[pid]]

# external force
set f_x    1.0
set f_y    0.0
set f_z    0.0

cellsystem domain_decomposition -no_verlet_list

# Constraints
set wall_left_id 1
set wall_right_id 2

constraint plane cell -10.0 -10.0 0.0 type $wall_left_id
constraint plane cell -10.0 -10.0 10.0 type $wall_right_id

#DPD-Thermostat
thermostat dpd $dpd_temperature $dpd_gamma $dpd_r_cut

# Solvent particles
for {set i 0} { $i < $n_solvent} { incr i } {
# Placing particles inside the constraints
    set posx [expr 1.0+8.0*[t_random]]
    set posy [expr 1.0+8.0*[t_random]]
    set posz [expr 1.0+8.0*[t_random]]

    set vx   [expr 0.1*[t_random]]
    set vy   [expr 0.1*[t_random]]
    set vz   [expr 0.1*[t_random]]
    
    part $i pos $posx $posy $posz type $solvent_id v $vx $vy $vz ext_force $f_x $f_y $f_z
}
galileiTransformParticles

# Interactions 
inter $wall_right_id $solvent_id lennard-jones $lj_eps $lj_sig $lj_cut $lj_shift $lj_off 
inter $wall_left_id $solvent_id lennard-jones $lj_eps $lj_sig $lj_cut $lj_shift $lj_off 
inter $wall_left_id $solvent_id tunable_slip $temp $gamma_L $r_cut_L $timestep $vx $vy $vz
inter $wall_right_id $solvent_id tunable_slip $temp $gamma_L $r_cut_L $timestep $vx $vy $vz

############################################
#      Procedures                          #
############################################

proc measure_kinetic_energy {} {
    global n_solvent
    global energy
    global n_solvent 
    global energy
    global E_ref

    set energy [analyze energy kinetic]
    set E_ref [expr $n_solvent*1.5]
    if {$energy <= $E_ref} {
	puts "Tunable-slip layer does not work ..."
	puts "Tunable slip boundaries fail!"
	exit 1;
    }
}

################### Integration #############################

for {set step 0} {$step < $int_loops} {incr step} {
    integrate $int_steps
}

measure_kinetic_energy

puts "Tunable-slip boundary conditions with constraints are ready..."
