// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef MMM1D_H
#define MMM1D_H

#include "config.h"
#include "particle_data.h"

#ifdef ELECTROSTATICS

typedef struct {
  double far_switch_radius_2;
  int    bessel_cutoff;
  int    bessel_calculated;
  double maxPWerror;
} MMM1D_struct;
extern MMM1D_struct mmm1d_params;

///
int printMMM1DToResult(Tcl_Interp *interp);

///
int inter_parse_mmm1d(Tcl_Interp *interp, int argc, char **argv);

///
int MMM1D_set_params(double switch_rad, int bessel_cutoff, double maxPWerror);

///
int MMM1D_tune(Tcl_Interp *interp);

///
void MMM1D_recalcTables();

///
int MMM1D_sanity_checks();

///
void MMM1D_init();

///
void calc_mmm1d_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist2,
				   double dist, double force[3]);

///
double mmm1d_coulomb_pair_energy(Particle *p1, Particle *p2, double d[3], double r2, double r);

#endif
#endif
