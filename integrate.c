/** \file integrate.c   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.h "integrate.h".
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "ghosts.h"
#include "debug.h"
#include "p3m.h"
#include "utils.h"
#include "thermostat.h"
#include "initialize.h"
#include "forces.h"

/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

double time_step = -1.0;
double sim_time = 0.0;

double skin = -1.0;
double max_range;
double max_range2;

int    particle_changed = 1;
int    interactions_changed = 1;
int    topology_changed = 1;
int    parameter_changed = 1;

double old_time_step;

double verlet_reuse=0.0;

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f]. */
void rescale_forces();
/** Propagate velocities: \f[ {\bf v} += {\bf f} \f]. */
void propagate_velocities();
/** Propagate positions: \f[ {\bf pos} += {\bf v} \f]. 
    Here also the verlet criterion is checked and communicated.
    \todo Put verlet criterion as inline function to verlet.h if possible.
*/
void propagate_positions(); 
/** combination of \ref propagate_velocities and \ref
    propagate_positions. */
void propagate_vel_pos();
/** combination of \ref rescale_forces and \ref
    propagate_velocities. */
void rescale_forces_propagate_vel(); 

/*@}*/

/************************************************************/

int invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  mpi_bcast_event(PARTICLE_CHANGED);
  mpi_bcast_event(INTERACTION_CHANGED);
  mpi_bcast_event(PARAMETER_CHANGED);
  mpi_bcast_event(TOPOLOGY_CHANGED);
  return TCL_OK;
}

int integrate(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <step num> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  /* translate argument */
  if (Tcl_GetInt(interp, argv[1], &n_steps) == TCL_ERROR)
    return (TCL_ERROR);

  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps", (char *) NULL);
    return (TCL_ERROR);
  }

  /* assume velocity verlet integration with langevin thermostat */
  if (argc != 2) {
    Tcl_AppendResult(interp, "too many arguments:  should be \"",
		     argv[0], " <task> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  mpi_integrate(n_steps);

  return (TCL_OK);
}

/************************************************************/

void integrate_vv_recalc_maxrange()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_recalc_maxrange:\n",this_node));

  /* sanity checks */
  if(time_step < 0.0 || skin < 0.0 || temperature < 0.0) {
    fprintf(stderr,"%d: ERROR: Can not initialize the integrator!:\n",this_node);
    if( time_step < 0.0 )
      fprintf(stderr,"%d: PROBLEM: You have to set the time_step!\n",this_node);
    if( skin < 0.0 )
      fprintf(stderr,"%d: PROBLEM: You have to set the skin!\n",this_node);
    if( temperature < 0.0 )
      fprintf(stderr,"%d: PROBLEM: You have to set the temperature!\n",this_node);
    errexit();
  }

  /* maximal interaction cutoff */
  calc_maximal_cutoff();
  max_range  = max_cut + skin;
  max_range2 = max_range * max_range;

  /* check real space interaction cutoff/range*/
  if(max_cut < 0.0) {
    fprintf(stderr,"%d: ERROR: You have to specify at least one interaction\n",this_node);
    errexit();
  }
  if((min_box_l/2.0) < max_range) {
    fprintf(stderr,"%d: ERROR: Maximal real space interaction %f is larger than half of the minimal box dimension %f\n",this_node,max_range,min_box_l);
    errexit();
  }

  if(min_local_box_l < max_range) {
    fprintf(stderr,"%d: ERROR: Maximal real space interaction %f is larger than minimal local box length %f\n",this_node,max_range,min_local_box_l);
    errexit();
  }
}

/************************************************************/

void integrate_vv(int n_steps)
{
  int i;
  int n_verlet_updates = 0;

  on_integration_start();
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps\n",this_node,
		      n_steps));

  if(parameter_changed || particle_changed || topology_changed || interactions_changed) {
    exchange_ghost();
    build_verlet_lists_and_force_calc();

    collect_ghost_forces();
    rescale_forces();
  }
    
  /* integration loop */
  INTEG_TRACE(fprintf(stderr,"%d START INTEGRATION: %d steps\n",this_node,n_steps));
  for(i=0;i<n_steps;i++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));
    propagate_vel_pos();
    if(rebuild_verletlist == 1) {
      INTEG_TRACE(fprintf(stderr,"%d: Rebuild Verlet List\n",this_node));
      n_verlet_updates++;
      invalidate_ghosts();
      exchange_and_sort_part();
      exchange_ghost();
      build_verlet_lists_and_force_calc();
    }
    else {
      update_ghost_pos();
      force_calc();
    }
    collect_ghost_forces();
    rescale_forces_propagate_vel();
    if(this_node==0) sim_time += time_step;
  }

  if(n_verlet_updates>0) verlet_reuse = n_steps/(double) n_verlet_updates;
  else verlet_reuse = 0;
  particle_changed     = 0; 
  interactions_changed = 0;
  topology_changed     = 0;
  parameter_changed    = 0;
}

/************************************************************/

void rescale_velocities() 
{
  Particle *p;
  int m,n,o,i, np;
  double scale;

  scale = time_step / old_time_step;
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].v[0] *= scale;
      p[i].v[1] *= scale;
      p[i].v[2] *= scale;
    }
  }
}

/* Callback functions */
/************************************************************/


int skin_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "skin must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  skin = data;
  mpi_bcast_parameter(FIELD_SKIN);
  mpi_bcast_event(PARAMETER_CHANGED);
  return (TCL_OK);
}

int time_step_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0.0) {
    Tcl_AppendResult(interp, "time step must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  old_time_step = time_step;
  time_step = data;
  mpi_set_time_step();

  mpi_bcast_event(PARAMETER_CHANGED);
  return (TCL_OK);
}

int time_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  sim_time = data;
  mpi_bcast_parameter(FIELD_SIMTIME);
  return (TCL_OK);
}


/* Privat functions */
/************************************************************/

void rescale_forces()
{
  Particle *p;
  int m,n,o,i, np;
  double scale;
  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].f[0] *= scale;
      p[i].f[1] *= scale;
      p[i].f[2] *= scale;

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f[0],p[i].f[1],p[i].f[2],p[i].v[0],p[i].v[1],p[i].v[2]));

    }
  }
}

void propagate_velocities() 
{
  Particle *p;
  int m,n,o,i, np;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_velocities:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{
	  p[i].v[0] += p[i].f[0];
	  p[i].v[1] += p[i].f[1];
	  p[i].v[2] += p[i].f[2];
	}
    }
  }
}

void rescale_forces_propagate_vel() 
{
  Particle *p;
  int m,n,o,i, np;
  double scale;

  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces_propagate_vel:\n",this_node));

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].f[0] *= scale;
      p[i].f[1] *= scale;
      p[i].f[2] *= scale;

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f[0],p[i].f[1],p[i].f[2],p[i].v[0],p[i].v[1],p[i].v[2]));

#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{
	  p[i].v[0] += p[i].f[0];
	  p[i].v[1] += p[i].f[1];
	  p[i].v[2] += p[i].f[2];
	}

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].v[0],p[i].v[1],p[i].v[2]));
      
    }
  }
}

void propagate_positions() 
{
  Particle *p;
  int m,n,o,i, np;

  int *verlet_flags = malloc(sizeof(int)*n_nodes);
  double skin2;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_positions:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin/2.0);

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{
	  p[i].r.p[0] += p[i].v[0];
	  p[i].r.p[1] += p[i].v[1];
	  p[i].r.p[2] += p[i].v[2];
	}
      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].p_old) > skin2 )
	rebuild_verletlist = 1; 
    }
  }

  /* communicate verlet criterion */
  MPI_Gather(&rebuild_verletlist, 1, MPI_INT, verlet_flags, 1, 
	     MPI_INT, 0, MPI_COMM_WORLD);
  if(this_node == 0)
    {
      rebuild_verletlist = 0;
      for(i=0;i<n_nodes;i++)
	if(verlet_flags[i]>0)
	  {
	    rebuild_verletlist = 1;
	    break;
	  }
    }
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Bcast(&rebuild_verletlist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  INTEG_TRACE(fprintf(stderr,"%d: prop_pos: rebuild_verletlist=%d\n",this_node,rebuild_verletlist));

  free(verlet_flags);
}

void propagate_vel_pos() 
{
  Particle *p;
  int m,n,o,i, np;
  int *verlet_flags = malloc(sizeof(int)*n_nodes);
  double skin2;

#ifdef ADDITIONAL_CHECKS
  double db_force,db_vel;
  double db_max_force=0.0, db_max_vel=0.0;
  int db_maxf_id=0,db_maxv_id=0;
#endif

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin/2.0);

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{
	  p[i].v[0] += p[i].f[0];
	  p[i].v[1] += p[i].f[1];
	  p[i].v[2] += p[i].f[2];
	

	  ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].v[0],p[i].v[1],p[i].v[2]));
	  
	  p[i].r.p[0] += p[i].v[0];
	  p[i].r.p[1] += p[i].v[1];
	  p[i].r.p[2] += p[i].v[2];
	}

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));


#ifdef ADDITIONAL_CHECKS
      /* force check */
      db_force = SQR(p[i].f[0])+SQR(p[i].f[1])+SQR(p[i].f[2]);
      if(db_force > skin2) 
	fprintf(stderr,"%d: Part %d has force %f (%f,%f,%f)\n",
		this_node,p[i].r.identity,sqrt(db_force),
		p[i].f[0],p[i].f[1],p[i].f[2]);
      if(db_force > db_max_force) { db_max_force=db_force; db_maxf_id=p[i].r.identity; }
      /* velocity check */
      db_vel   = SQR(p[i].v[0])+SQR(p[i].v[1])+SQR(p[i].v[2]);
      if(db_vel > skin2) 
	fprintf(stderr,"%d: Part %d has velocity %f (%f,%f,%f)\n",
		this_node,p[i].r.identity,sqrt(db_vel),
		p[i].v[0],p[i].v[1],p[i].v[2]);
      if(db_vel > db_max_vel) { db_max_vel=db_vel; db_maxv_id=p[i].r.identity; }
#endif

      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].p_old) > skin2 )
	rebuild_verletlist = 1; 
    }
  }

#ifdef ADDITIONAL_CHECKS
  if(db_max_force > skin2) 
    fprintf(stderr,"%d: max_force=%e, part=%d f=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_force),db_maxf_id,local_particles[db_maxf_id]->f[0],
	    local_particles[db_maxf_id]->f[1],local_particles[db_maxf_id]->f[2]);
  if(db_max_vel > skin2)
    fprintf(stderr,"%d: max_vel=%e, part=%d v=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_vel),db_maxv_id,local_particles[db_maxv_id]->v[0],
	    local_particles[db_maxv_id]->v[1],local_particles[db_maxv_id]->v[2]);
#endif

  /* communicate verlet criterion */
  MPI_Gather(&rebuild_verletlist, 1, MPI_INT, verlet_flags, 1, 
	     MPI_INT, 0, MPI_COMM_WORLD);
  if(this_node == 0)
    {
      rebuild_verletlist = 0;
      for(i=0;i<n_nodes;i++)
	if(verlet_flags[i]>0)
	  {
	    rebuild_verletlist = 1;
	    break;
	  }
    }
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Bcast(&rebuild_verletlist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  INTEG_TRACE(fprintf(stderr,"%d: prop_pos: rebuild_verletlist=%d\n",this_node,rebuild_verletlist));

  free(verlet_flags);
}

