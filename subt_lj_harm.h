// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef SUBT_LJ_HARM_H
#define SUBT_LJ_HARM_H
/** \file subt_lj_harm.h
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force 
 *  from the HARMONIC Energy and/or the force for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sathish@mpip-mainz.mpg.de">sathish</a>
*/

#ifdef LENNARD_JONES

/************************************************************/

/// parameters for the subtract lj from a harmonic bond potential
MDINLINE int subt_lj_harm_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.subt_lj_harm.k = k;
  bonded_ia_params[bond_type].p.subt_lj_harm.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_SUBT_LJ_HARM;
  bonded_ia_params[bond_type].p.subt_lj_harm.r2 = SQR(bonded_ia_params[bond_type].p.subt_lj_harm.r);
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/** Computes the difference between the HARMONIC and the LENNARD-JONES pair forces 
    and adds this force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of this interaction (see \ref #inter).
*/
MDINLINE int calc_subt_lj_harm_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac_harm=0.0, fac_lj=0.0, fac;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  fac_harm = iaparams->p.subt_lj_harm.k;
  fac_harm *= (dist-iaparams->p.subt_lj_harm.r);
  fac_harm /= dist;

  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;

    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);			  
    }
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
    }      
  }    

  fac = -(fac_harm + fac_lj);
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_HARM f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_HARM f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

  return 0;
}

MDINLINE int subt_lj_harm_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double energy_harm=0.0, energy_lj=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  energy_harm = 0.5*iaparams->p.subt_lj_harm.k;
  energy_harm *= SQR(dist-iaparams->p.subt_lj_harm.r);
  
  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
    r_off = dist - ia_params->LJ_offset;
    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }   
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  *_energy = energy_harm-energy_lj;

  return 0;
}

#endif
#endif
