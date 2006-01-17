// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file blockfile_tcl.h
    Contains only the tcl interface for block coded files.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    It is the header file for \ref blockfile_tcl.c "blockfile_tcl.c" and provides the
    tcl command \ref tcl_blockfile.
*/
#ifndef BLOCKFILE_TCL_H
#define BLOCKFILE_TCL_H
#include <tcl.h>

/** Implementation of the Tcl command \ref tcl_blockfile. Allows to read and write
    blockfile comfortably from Tcl. */
int blockfile(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
#endif
