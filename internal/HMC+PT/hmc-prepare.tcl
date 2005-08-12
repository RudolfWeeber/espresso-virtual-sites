#!/bin/sh
#This script is used to create the initial starting configurations for all ensembles.

source ./hmc-initialize.tcl
source ./hmc-setup.tcl

for {set ens_i 0} { $ens_i < $ensemble_num } {incr ens_i} {
    set bjerrum [lindex $bjerrum_list $ens_i]
    exec mkdir  lb_[format %06.3f $bjerrum]
    exec cp "equilibrized.gz"  "lb_[format %06.3f $bjerrum]/HMC.t_00000.gz"
    set lockf [open "lb_[format %06.3f $bjerrum]/swpd_[format %05d -1]" "w"]; puts $lockf " "; close $lockf
}

set f [open "wherewasI" "w"]
puts $f "0"
close $f
#exit 0
#############################################################
#      Set up the coulomb information
#############################################################
#use the maximum value of bjerrum to setup the coulomb interactions
puts "Setting up the coulomb interaction"
set coulombfile [open "coulomb.txt" "w"]
puts $coulombfile "[inter coulomb $bjerrum p3m tune accuracy 0.0001]"
close $coulombfile


exit 0
