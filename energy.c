// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file energy.c
    Implementation of \ref energy.h "energy.h".
*/
#include "energy.h"
#include "parser.h"
#include "cells.h"
#include "integrate.h"
#include "initialize.h"
#include "domain_decomposition.h"
#include "nsquare.h"
#include "layered.h"
#include "elc.h"

Observable_stat energy = {0, {NULL,0,0}, 0,0,0};
Observable_stat total_energy = {0, {NULL,0,0}, 0,0,0};

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range energies (P3M, EWALD, MMM2d...). */
void calc_long_range_energies();

/** allocate energy arrays and initialize with zero */
void init_energies(Observable_stat *stat);

/** on the master node: calc energies only if necessary */
void master_energy_calc();

/************************************************************/

void energy_calc(double *result)
{
  if (!check_obs_calc_initialized())
    return;

  init_energies(&energy);

  on_observable_calc();
  
  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_energies();
    break;
  case CELL_STRUCTURE_DOMDEC: 
    if(dd.use_vList) {
      if (rebuild_verletlist)  
	build_verlet_lists();
      calculate_verlet_energies();
    }
    else
      calculate_link_cell_energies();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_energies();
  }
  /* rescale kinetic energy */
  energy.data.e[0] /= (2.0*time_step*time_step);

  calc_long_range_energies();
  
  /* gather data */
  MPI_Reduce(energy.data.e, result, energy.data.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/************************************************************/

void calc_long_range_energies()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
  case COULOMB_P3M:
    P3M_charge_assign(); 
    energy.coulomb[1] = P3M_calc_kspace_forces(0,1);
    break;
  case COULOMB_ELC_P3M:
    // assign the original charges first
    // they may not have been assigned yet
    P3M_charge_assign(); 
    if(!elc_params.dielectric_contrast_on)
      energy.coulomb[1] = P3M_calc_kspace_forces(0,1);
    else {
      energy.coulomb[1] = 0.5*P3M_calc_kspace_forces(0,1); 
      energy.coulomb[1]+= 0.5*ELC_P3M_dielectric_layers_energy_self(); 

      //  assign both original and image charges now
      ELC_P3M_charge_assign_both();
      ELC_P3M_modify_p3m_sums_both();

      energy.coulomb[1] += 0.5*P3M_calc_kspace_forces(0,1); 

      //assign only the image charges now
      ELC_P3M_charge_assign_image();
      ELC_P3M_modify_p3m_sums_image();

      energy.coulomb[1]-= 0.5*P3M_calc_kspace_forces(0,1); 
    }
    energy.coulomb[2] = ELC_energy();
    break;
  case COULOMB_EWALD:
    energy.coulomb[1] = EWALD_calc_kspace_forces(0,1);
    EWALD_TRACE(fprintf(stderr,"%d: EWALD: energy.coulomb[1]=%g\n",this_node,energy.coulomb[1]));
    break;
  case COULOMB_MMM2D:
    *energy.coulomb += MMM2D_far_energy();
    *energy.coulomb += MMM2D_dielectric_layers_energy_contribution();
    break;
    /* calculate electric part of energy (only for MAGGS) */
  case COULOMB_MAGGS:
    *energy.coulomb += maggs_electric_energy();
    break;
  }
#endif
}

/************************************************************/

void init_energies(Observable_stat *stat)
{
  int n_pre, n_non_bonded, n_coulomb;

  n_pre        = 1;
  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  n_coulomb    = 0;
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:  n_coulomb = 0; break;
  case COULOMB_ELC_P3M: n_coulomb = 3; break;
  case COULOMB_P3M:   n_coulomb = 2; break;
  case COULOMB_EWALD: n_coulomb = 2; break;
  default: n_coulomb  = 1;
  }
#endif
  
  obsstat_realloc_and_clear(stat, n_pre, n_bonded_ia, n_non_bonded, n_coulomb, 1);
  stat->init_status = 0;
}

/************************************************************/

void master_energy_calc() {
  mpi_gather_stats(1, total_energy.data.e, NULL, NULL, NULL);

  total_energy.init_status=1;
}

/****************************************************************************************
 *                                 parser
 ****************************************************************************************/

static void print_detailed_energies(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value;
  int i, j;

  value = total_energy.data.e[0];
  for (i = 1; i < total_energy.data.n; i++)
    value += total_energy.data.e[i];

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, "{ energy ", buffer, " } ", (char *)NULL);

  Tcl_PrintDouble(interp, total_energy.data.e[0], buffer);
  Tcl_AppendResult(interp, "{ kinetic ", buffer, " } ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
      Tcl_PrintDouble(interp, *obsstat_bonded(&total_energy, i), buffer);
      Tcl_AppendResult(interp,
		       get_name_of_bonded_ia(bonded_ia_params[i].type),
		       " ", buffer, " } ", (char *) NULL);
    }
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded(&total_energy, i, j), buffer);
	Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);	    
      }
    }

#ifdef ELECTROSTATICS
  if(coulomb.method != COULOMB_NONE) {
    /* total Coulomb energy */
    value = total_energy.coulomb[0];
    for (i = 1; i < total_energy.n_coulomb; i++)
      value += total_energy.coulomb[i];
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);

    /* if it is split up, then print the split up parts */
    if (total_energy.n_coulomb > 1) {
      for (i = 0; i < total_energy.n_coulomb; i++) {
	Tcl_PrintDouble(interp, total_energy.coulomb[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
    Tcl_AppendResult(interp, " }", (char *)NULL);
  }
#endif
}

/************************************************************/

int parse_and_print_energy(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze energy [{ fene <type_num> | harmonic <type_num> | subt_lj_harm <type_num> | subt_lj_fene <type_num> | subt_lj <type_num> | lj <type1> <type2> | ljcos <type1> <type2> | ljcos2 <type1> <type2> | gb <type1> <type2> | coulomb | kinetic | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j;
  double value;
  value = 0.0;
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  if (total_energy.init_status == 0) {
    init_energies(&total_energy);
    master_energy_calc();
  }

  if (argc == 0)
    print_detailed_energies(interp);
  else {

    if      (ARG0_IS_S("kinetic"))
      value = total_energy.data.e[0];
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze energy bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_bonded(&total_energy, i);
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
	     ARG0_IS_S("lj-cos") ||
             ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("gb") ||
	     ARG0_IS_S("tabulated")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze energy nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded(&total_energy, i, j);
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      value = total_energy.coulomb[0];
      for (i = 1; i < total_energy.n_coulomb; i++)
	value += total_energy.coulomb[i];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if (ARG0_IS_S("total")) {
      value = total_energy.data.e[0];
      for (i = 1; i < total_energy.data.n; i++)
	value += total_energy.data.e[i];
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze energy",
		       (char *)NULL);
      return (TCL_ERROR);
    }
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }

  return (TCL_OK);
}
