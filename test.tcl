#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; export EF_ALLOW_MALLOC_0=1; if [ "$1" != "" ]; then NP=$1; else NP=8; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else lamboot -v hostfile; exec mpirun -c2c -toff -O -np $NP -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

##################################################
# settings
##################################################
set tcl_precision 5

set intsteps 100
set maxtime  50

setmd periodic 1 1 1
setmd bjerrum 1.0
setmd p3m_alpha 0.27
setmd p3m_r_cut 3.0
setmd p3m_mesh 8 8 8
setmd p3m_cao 3 5000
setmd p3m_epsilon 0.1
setmd p3m_mesh_off 0.5 0.5 0.5

setmd box_l 10.0 10.0 10.0
setmd max_num_cells 512
# integrator
setmd skin 0.4
setmd time_step 0.0001

# read particles
##################################################

if { ![file exists "config.gz"]} {
    error "please generate a configuration file with setup.tcl"
    exit
}
set f [open "|gzip -cd config.gz" r]
while {[blockfile $f read auto] != "eof" } {}
close $f
puts "n_part = [part number]"
puts "grid   = \{[setmd node_grid]\}"
puts "box    = \{[setmd box]\}"

# setup interactions
##################################################
inter 0 0 lennard-jones 1 1 1.12246 0 0
inter 1 1 lennard-jones 1 1 1.12246 0 0
inter 2 2 lennard-jones 1 1 1.12246 0 0
inter 0 1 lennard-jones 1 1 1.12246 0 0
inter 0 2 lennard-jones 1 1 1.12246 0 0
inter 1 2 lennard-jones 1 1 1.12246 0 0

inter 1 fene 1 2

puts "add first bond"
part 1 bond 1 3
puts [part 1]
part 1 bond delete 1 3
puts [part 1]

# integration
##################################################

#for {set port 10000} { $port < 65000 } { incr port } {
#    catch {imd connect $port} res
#    if {$res == ""} break
#}
#puts "opened port $port"

integrate init

puts [time {
  for {set i 0} { $i < $maxtime } { incr i } {
      puts "step $i"
      integrate $intsteps
  #   imd pos
  puts "mindist=[mindist]"
  }
}]

integrate exit

#puts "imd: [imd disconnect]"

# exit
##################################################
puts "finished"
exit
