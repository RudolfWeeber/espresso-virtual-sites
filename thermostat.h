// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef THERMOSTAT_H
#define THERMOSTAT_H
/** \file thermostat.h 


    <b>Responsible:</b>
    <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>

*/

#include <tcl.h>
#include <math.h>
#include "utils.h"
#include "particle_data.h"
#include "parser.h"
#include "random.h"
#include "global.h"
#include "communication.h"
#include "integrate.h"
#include "cells.h"
#include "pressure.h"
#include "lb.h"

/** \name Thermostat switches*/
/************************************************************/
/*@{*/

#define THERMO_OFF      0
#define THERMO_LANGEVIN 1
#define THERMO_DPD      2
#define THERMO_NPT_ISO  4
#define THERMO_LB       8

/*@}*/

/************************************************
 * exported variables
 ************************************************/

/** Switch determining which thermostat to use. This is a or'd value
    of the different possible thermostats (defines: \ref THERMO_OFF,
    \ref THERMO_LANGEVIN, \ref THERMO_DPD \ref THERMO_NPT_ISO). If it
    is zero all thermostats are switched off and the temperature is
    set to zero.  */
extern int thermo_switch;

/** temperature. */
extern double temperature;

/** Langevin friction coefficient gamma. */
extern double langevin_gamma;

/** DPD Friction coefficient gamma. */
extern double dpd_gamma;
/** DPD thermostat cutoff */
extern double dpd_r_cut;
/** DPD thermostat weight function */
extern int dpd_wf;
#ifdef TRANS_DPD
/** DPD transversal Friction coefficient gamma. */
extern double dpd_tgamma;
/** trans DPD thermostat cutoff */
extern double dpd_tr_cut;
/** trans DPD thermostat weight function */
extern int dpd_twf;
#endif

/** Friction coefficient for nptiso-thermostat's inline-function friction_therm0_nptiso */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's inline-function friction_thermV_nptiso */
extern double nptiso_gammav;

/************************************************
 * functions
 ************************************************/

/** Implementation of the tcl command \ref tcl_thermostat. This function
    allows to change the used thermostat and to set its parameters.
 */
int thermostat(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** initialize constants of the thermostat on
    start of integration */
void thermo_init();

/** very nasty: if we recalculate force when leaving/reentering the integrator,
    a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
    numbers are drawn twice, resulting in a different variance of the random force.
    This is corrected by additional heat when restarting the integrator here.
    Currently only works for the Langevin thermostat, although probably also others
    are affected.
*/
void thermo_heat_up();

/** pendant to \ref thermo_heat_up */
void thermo_cool_down();

#ifdef NPT
/** add velocity-dependend noise and friction for NpT-sims to the particle's velocity 
    @param dt_vj  j-component of the velocity scaled by time_step dt 
    @return       j-component of the noise added to the velocity, also scaled by dt (contained in prefactors) */
MDINLINE double friction_therm0_nptiso(double dt_vj) {
  extern double nptiso_pref1, nptiso_pref2;
  if(thermo_switch & THERMO_NPT_ISO)   
    return ( nptiso_pref1*dt_vj + nptiso_pref2*(d_random()-0.5) );
  return 0.0;
}

/** add p_diff-dependend noise and friction for NpT-sims to \ref nptiso_struct::p_diff */
MDINLINE double friction_thermV_nptiso(double p_diff) {
  extern double nptiso_pref3, nptiso_pref4;
  if(thermo_switch & THERMO_NPT_ISO)   
    return ( nptiso_pref3*p_diff + nptiso_pref4*(d_random()-0.5) );
  return 0.0;
}
#endif

/** Callback marking setting the temperature as outdated */
int thermo_ro_callback(Tcl_Interp *interp, void *_data);

/** overwrite the forces of a particle with
    the friction term, i.e. \f$ F_i= -\gamma v_i + \xi_i\f$.
*/
MDINLINE void friction_thermo_langevin(Particle *p)
{
  extern double langevin_pref1, langevin_pref2;

  int j;
#ifdef MASS
  double massf = sqrt(PMASS(*p));
#else
  double massf = 1;
#endif

  for ( j = 0 ; j < 3 ; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p->l.ext_flag & COORD_FIXED(j)))
#endif
      p->f.f[j] = langevin_pref1*p->m.v[j]*PMASS(*p) + langevin_pref2*(d_random()-0.5)*massf;
#ifdef EXTERNAL_FORCES
    else p->f.f[j] = 0;
#endif
  }
  

  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
  THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
}

#ifdef ROTATION
/** set the particle torques to the friction term, i.e. \f$\tau_i=-\gamma w_i + \xi_i\f$.
    The same friction coefficient \f$\gamma\f$ is used as that for translation.
*/
MDINLINE void friction_thermo_langevin_rotation(Particle *p)
{
  extern double langevin_pref2;

  int j;
#ifdef EXTERNAL_FORCES
  if(!(p->l.ext_flag & COORDS_FIX_MASK))
#endif
    {
      for ( j = 0 ; j < 3 ; j++)
	p->f.torque[j] = -langevin_gamma*p->m.omega[j] + langevin_pref2*(d_random()-0.5);

      ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
      THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
    }
}
#endif

#ifdef DPD
/** Calculate Random Force and Friction Force acting between particle
    p1 and p2 and add them to their forces. */
MDINLINE void add_dpd_thermo_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double dist2)
{
  extern double dpd_gamma,dpd_pref1, dpd_pref2,dpd_r_cut,dpd_r_cut_inv;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma, dpd_pref3, dpd_pref4,dpd_tr_cut,dpd_tr_cut_inv;
  extern int dpd_twf;
#endif
  int j;
  // velocity difference between p1 and p2
  double vel12_dot_d12=0.0;
  // inverse distance
  double dist_inv;
  // weighting functions for friction and random force
  double omega,omega2;// omega = w_R/dist
  double friction, noise;
  //Projection martix
#ifdef TRANS_DPD
  int i;
  double P_times_dist_sqr[3][3]={{dist2,0,0},{0,dist2,0},{0,0,dist2}},noise_vec[3];
  double f_D[3],f_R[3];
#endif
  double tmp;
  dist_inv = 1.0/dist;
  if((dist < dpd_r_cut)&&(dpd_gamma > 0.0)) {
    if ( dpd_wf == 1 )
    {
       omega    = dist_inv;
    }
    else 
    {
    	omega    = dist_inv- dpd_r_cut_inv;
    }
    omega2   = SQR(omega);
    //DPD part
    if (dpd_gamma > 0.0 ){
       // friction force prefactor
      for(j=0; j<3; j++)  vel12_dot_d12 += (p1->m.v[j] - p2->m.v[j]) * d[j];
      friction = dpd_pref1 * omega2 * vel12_dot_d12;
      // random force prefactor
      noise    = dpd_pref2 * omega      * (d_random()-0.5);
      for(j=0; j<3; j++) {
        p1->f.f[j] += ( tmp = (noise - friction)*d[j] );
        p2->f.f[j] -= tmp;
      }
    }
  }
#ifdef TRANS_DPD
    //DPD2 part
  if ((dist < dpd_tr_cut)&&(dpd_tgamma > 0.0)){
      if ( dpd_twf == 1 )
      {
        omega    = dist_inv;
      }
      else 
      {
        omega    = dist_inv- dpd_tr_cut_inv;
      }
      omega2   = SQR(omega);
      for (i=0;i<3;i++){
        //noise vector
        noise_vec[i]=d_random()-0.5;
        // Projection Matrix
        for (j=0;j<3;j++){
          P_times_dist_sqr[i][j]-=d[i]*d[j];
        }
      }
      for (i=0;i<3;i++){
        //Damping force
        f_D[i]=0;
        //Random force
        f_R[i]=0;
        for (j=0;j<3;j++){
          f_D[i]+=P_times_dist_sqr[i][j]*(p1->m.v[j] - p2->m.v[j]);
          f_R[i]+=P_times_dist_sqr[i][j]*noise_vec[j];
        }
        f_D[i]*=dpd_pref3*omega2;
        f_R[i]*=dpd_pref4*omega*dist_inv;
      }
      for(j=0; j<3; j++) {
        tmp=f_R[j]-f_D[j];
        p1->f.f[j] += tmp;
        p2->f.f[j] -= tmp;
      }
  }
#endif
}
#endif

#ifdef LB
int thermo_parse_lb(Tcl_Interp * interp, int argc, char ** argv);
#endif

#endif
