#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=1; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Demo of a buckyball                                      #
#                                                           #
#                                                           #
#  Created:       26.08.2003 by AxA & BAM                   #
#                                                           #
#############################################################

set case2timeout 120

#############################################################
#  Setup GUI                                                #
#############################################################

# popup window if desired
proc popup {text} {
    toplevel .tag
    button .tag.btn -text "text" -command "destroy .tag"
    pack .tag.btn
    wm geometry .tag +60+0
    tkwait window .tag
    update
}

##### case 1
frame .case1 -relief raised -border 2
pack .case1 -expand 1 -fill both -in .

# title
frame .case1.title
label .case1.title.l -justify left -text "Schrumpfen"
button .case1.title.b -text "Start" -command Case1Start
pack .case1.title.l -fill both -side left -in .case1.title
pack .case1.title.b -side right -in .case1.title
pack .case1.title -pady {0 10} -fill x -side top -in .case1

# sliders
label .case1.label -justify left -text "St�rke der Elektrostatik"
scale .case1.slider -orient h -from 0 -to 15 \
    -resolution [expr 10/100.] -command Case1BjerrumChange
pack .case1.slider .case1.label -fill x -in .case1

label .case1.label2 -justify left -text "Temperatur"
scale .case1.slider2 -orient h -from 0 -to 2 \
    -resolution [expr 10/100.] -command Case1TempChange
pack .case1.slider2 .case1.label2 -fill x -in .case1

set disabledfg [.case1.label cget -disabledforeground]
set normalfg  [.case1.label cget -foreground]

##### case 2
frame .case2 -relief raised -border 2
pack .case2 -expand 1 -fill both -in .

# title
frame .case2.title
label .case2.title.l -justify left -text "Irrgarten"
button .case2.title.b -text "Start" -command Case2Start
pack .case2.title.l -fill both -side left -in .case2.title
pack .case2.title.b -side right -in .case2.title
pack .case2.title -pady {0 10} -fill x -side top -in .case2

# sliders
label .case2.label1 -justify left -text "Ladung Wand 1"
scale .case2.slider1 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall1ChargeChange
pack .case2.slider1 .case2.label1 -fill x -in .case2

label .case2.label2 -justify left -text "Ladung Wand 2"
scale .case2.slider2 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall2ChargeChange
pack .case2.slider2 .case2.label2 -fill x -in .case2

label .case2.label3 -justify left -text "Ladung Wand 3"
scale .case2.slider3 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall3ChargeChange
pack .case2.slider3 .case2.label3 -fill x -in .case2

label .case2.status -justify left -text ""
pack .case2.status -fill x -in .case2

proc disableCase1 {} {
    global disabledfg
    .case1.label configure -state disabled
    .case1.slider configure -state disabled -foreground $disabledfg 
    .case1.label2 configure -state disabled
    .case1.slider2 configure -state disabled -foreground $disabledfg 
}

proc enableCase1 {} {
    global normalfg
    .case1.label configure -state normal
    .case1.slider configure -state normal -foreground $normalfg 
    .case1.slider set 8
    .case1.label2 configure -state normal
    .case1.slider2 configure -state normal -foreground $normalfg 
    .case1.slider2 set 0

    .case2.status configure -text ""
}

proc disableCase2 {} {
    global disabledfg displayclock
    .case2.label1 configure -state disabled
    .case2.slider1 configure -state disabled -foreground $disabledfg
    .case2.label2 configure -state disabled
    .case2.slider2 configure -state disabled -foreground $disabledfg
    .case2.label3 configure -state disabled
    .case2.slider3 configure -state disabled -foreground $disabledfg

    .case2.status configure -text ""
    set displayclock 0
}

proc enableCase2 {} {
    global normalfg
    .case2.label1 configure -state normal
    .case2.slider1 configure -state normal -foreground $normalfg
    .case2.slider1 set 0
    .case2.label2 configure -state normal
    .case2.slider2 configure -state normal -foreground $normalfg
    .case2.slider2 set 0
    .case2.label3 configure -state normal
    .case2.slider3 configure -state normal -foreground $normalfg
    .case2.slider3 set 0
}

proc disableStarts {} {
    .case1.title.b configure -state disabled
    .case2.title.b configure -state disabled
}

proc enableStarts {} {
    .case1.title.b configure -state normal
    .case2.title.b configure -state normal
}

#############################################################
#      Helpers                                              #
#############################################################

proc imd_reconnect {case} {
#    global vmd_pid

    imd disconnect

    catch {exec killall vmd_LINUX}

    while { [catch {imd connect 10000} res] } {
	.case2.status configure -text "Waiting for port to unbind...\nerror: $res"
	update
	after 500
    }
    after 1000
    exec vmd -e $case.script &
}

#############################################################
#      Events                                               #
#############################################################

proc Case1Start {} {
    global run_sim

    .case2.status configure -text "Starte..."    
    set run_sim 0
    update

    disableStarts
    disableCase1
    disableCase2

    part deleteall
    set inp [open "case1.blk" r]
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp
    puts "Read [setmd n_part] particles..."

    imd_reconnect case1

    after 3000

    enableCase1
    enableStarts
    set run_sim 1

    .case2.status configure -text "Fertig"
}

set displayclock 0
proc Case2Start {} {
    global run_sim case2stime displayclock N0 N1 N2 N3

    .case2.status configure -text "Starte..."    
    set run_sim 0
    update

    disableStarts
    disableCase1
    disableCase2

    part deleteall
    set inp [open "case2.blk" r]
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp
    for {set i $N0} {$i < [expr $N0+$N1+$N2+$N3]} {incr i} { part $i fix }
    puts "Read & fixed [setmd n_part] particles..."

    imd_reconnect case2

    after 3000
    .case2.status configure -text "...3..."
    update; after 1000
    .case2.status configure -text "...2..."
    update; after 1000
    .case2.status configure -text "...1..."
    update; after 1000

    enableCase2
    enableStarts
    set run_sim 2

    .case2.status configure -text "Los!"
    set case2stime [clock seconds]
    set displayclock 1
}

proc Case1BjerrumChange {l} {
    global run_sim
    if { $run_sim == 1} { 
	inter coulomb $l dh 0 [expr 2*$l]
	puts "Changed Bjerrum length to $l."
    }
}

proc Case1TempChange {T} {
    global run_sim
    if { $run_sim == 1} { 
	setmd temp $T
	puts "Changed temperature to $T."
    }
}

proc Case2Wall1ChargeChange {c} {
    global run_sim N0 N1
    if { $run_sim == 2} { 
	for {set i [expr $N0]} {$i < [expr $N0+$N1]} {incr i} { part $i q $c }
	puts "Charged Wall 1 with $c."
    }
}

proc Case2Wall2ChargeChange {c} {
    global run_sim N0 N1 N2
    if { $run_sim == 2 } { 
	for {set i [expr $N0+$N1]} {$i < [expr $N0+$N1+$N2]} {incr i} { part $i q $c } 
	puts "Charged Wall 2 with $c."
    }
}

proc Case2Wall3ChargeChange {c} {
    global run_sim N0 N1 N2 N3
    if { $run_sim } { 
	for {set i [expr $N0+$N1+$N2]} {$i < [expr $N0+$N1+$N2+$N3]} {incr i} { part $i q $c }
	puts "Charged Wall 3 with $c."
    }
}

#############################################################
#      Integration                                          #
#############################################################

disableCase1
disableCase2

set displayclock 0
set run_sim 0

while { 1 } {
    if { $run_sim != 0 } {
#	set c1 [part 1 print pos]; set c2 [part 2 print pos]
#	set x  [expr [lindex $c1 0]-[lindex $c2 0]]
#	set y  [expr [lindex $c1 1]-[lindex $c2 1]]
#	set z  [expr [lindex $c1 2]-[lindex $c2 2]]
#	puts "bond_l between 1 and 2 = [expr sqrt($x*$x + $y*$y + $z*$z)]"
	integrate 100
	catch [imd positions]
    } {
	after 100
	imd listen 1
    }
    if { $displayclock } {
	set etime [expr [clock seconds] - $case2stime]
	.case2.status configure -text [clock format $etime -format "%M:%S"]
	if { $etime > $case2timeout } {
	    .case2.status configure -text "Verloren!"
	    set displayclock 0
	    set run_sim 0
	}
    }
    update
}