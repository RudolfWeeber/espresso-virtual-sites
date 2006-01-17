// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file mmm-common.h
    modified polygamma functions. See Arnold,Holm 2002
*/
#ifndef MMM_COMMON_H
#define MMM_COMMON_H

#include "polynom.h"
#include "specfunc.h"

/** \name Math Constants */
/*@{*/
#define C_2PI     (2*M_PI)
#define C_GAMMA   0.57721566490153286060651209008
#define C_2LOG4PI -5.0620484939385815859557831885
#define C_2PISQR  C_2PI*C_2PI
/*@}*/

/** table of the Taylor expansions of the modified polygamma functions */
extern Polynom *modPsi;
extern int      n_modPsi;

/** modified polygamma for even order 2*n, n >= 0 */
MDINLINE double mod_psi_even(int n, double x)
{ return evaluateAsTaylorSeriesAt(&modPsi[2*n],x*x); }

/** modified polygamma for odd order 2*n+1, n>= 0 */
MDINLINE double mod_psi_odd(int n, double x)
{ return x*evaluateAsTaylorSeriesAt(&modPsi[2*n+1], x*x); }

/** create the both the even and odd polygamma functions up to order 2*n */
void create_mod_psi_up_to(int n);

#endif
