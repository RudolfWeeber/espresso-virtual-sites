// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef CONSTRAINT_H
#define CONSTRAINT_H

/** \file constraint.h
 *  Routines for handling of constraints. Implemented are walls, cylinders and spheres.
 *  Only activ if you define CONSTRAINTS in \ref config.h.
 *  see also \ref interaction_data.h
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 */

#ifdef CONSTRAINTS

MDINLINE void calculate_wall_dist(Particle *p1, Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;
  IA_parameters *ia_params;

  ia_params=get_ia_param(p1->r.type, c_p->r.type);

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += p1->r.p[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
}


MDINLINE void calculate_sphere_dist(Particle *p1, Particle *c_p, Constraint_sphere *c, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;
  IA_parameters *ia_params;
 
  ia_params=get_ia_param(p1->r.type, c_p->r.type);
 
  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = c->pos[i] - p1->r.p[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = sqrt(c_dist);
  
  if ( c->direction == -1 ) {
  /* apply force towards inside the sphere */
    *dist = c->rad - c_dist;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= fac;}
  else  {
   /* apply force towards outside the sphere */
    *dist = c_dist - c->rad;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= -fac;
  }
}


MDINLINE void calculate_maze_dist(Particle *p1, Particle *c_p, Constraint_maze *c, double *dist, double *vec)
{
  int i,min_axis,cursph[3],dim;
  float diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
  float sph_vec[3],cyl_vec[3];
  IA_parameters *ia_params;
 
  ia_params = get_ia_param(p1->r.type, c_p->r.type);
  dim=(int) c->dim;
  diasph = box_l[0]/c->nsphere;
  
  /* First determine the distance to the sphere */
  c_dist=0.0;
  for(i=0;i<3;i++) {
    cursph[i] = (int) (p1->r.p[i]/diasph);
    sph_vec[i] = (cursph[i]+0.5) * diasph  - (p1->r.p[i]);
    c_dist += SQR(sph_vec[i]);
  }
  c_dist = sqrt(c_dist);
  sph_dist = c->sphrad - c_dist;
  fac = sph_dist / c_dist;
  for(i=0;i<3;i++) cyl_vec[i] = sph_vec[i];
  for(i=0;i<3;i++) sph_vec[i] *= fac;
  
  /* Now calculate the cylinder stuff */
  /* have to check all 3 cylinders */
  min_axis=2;
  cyl_dist=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[1]*cyl_vec[1]);
  
  if(dim > 0 ){
    temp_dis=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[2]*cyl_vec[2]);
    if ( temp_dis < cyl_dist) {
        min_axis=1;
        cyl_dist=temp_dis;
    }

    if(dim > 1 ){
        temp_dis=sqrt(cyl_vec[1]*cyl_vec[1]+cyl_vec[2]*cyl_vec[2]);
        if ( temp_dis < cyl_dist) {
            min_axis=0;
            cyl_dist=temp_dis;
        }
    }
  }
  cyl_vec[min_axis]=0.;
  
  c_dist=cyl_dist;
  cyl_dist = c->cylrad - c_dist;
  fac = cyl_dist / c_dist;
  for(i=0;i<3;i++) cyl_vec[i] *= fac;
  
  /* Now decide between cylinder and sphere */
  if ( sph_dist > 0 ) {
    if ( sph_dist>cyl_dist ) {
        *dist = sph_dist;
        for(i=0;i<3;i++) vec[i] = sph_vec[i];
    } else {
        *dist = cyl_dist;
        for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
    }
  } else {
    *dist = cyl_dist;
    for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
  }
}

MDINLINE void calculate_cylinder_dist(Particle *p1, Particle *c_p, Constraint_cylinder *c, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, c_p->r.type);


  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = p1->r.p[i] - c->pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * c->axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * c->axis[i] ;
    d_per_vec[i] = p1->r.p[i] - (c->pos[i] + d_par_vec[i]) ;
  }
		
  d_per=sqrt(SQR(d_real)-SQR(d_par));
  d_par = fabs(d_par) ;
    
  if ( c->direction == -1 ) {
    /*apply force towards inside cylinder */
    d_per = c->rad - d_per ;
    d_par = c->length - d_par;
    if (d_per < d_par )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= d_per_vec[i] * d_per /  (c->rad - d_per) ;
      }
    } else {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= d_par_vec[i] * d_par /  (c->length - d_par) ;
      }
    }
  } else {
    /*apply force towards outside cylinder */
    d_per = d_per - c->rad ;
    d_par = d_par - c->length ;
    if (d_par < 0 )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= - d_per_vec[i] * d_per /  (d_per + c->rad) ;
      }
    } else if ( d_per < 0) {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= - d_par_vec[i] * d_par /  (d_par + c->length) ;
      }
    } else {
      *dist = sqrt( SQR(d_par) + SQR(d_per)) ;
      for (i=0; i<3;i++) {
	vec[i]= - d_per_vec[i] * d_per /  (d_per + c->rad) ;
	vec[i]+= - d_par_vec[i] * d_par /  (d_par + c->length) ;
      }			
    }
  }
}

MDINLINE void add_rod_force(Particle *p1, Particle *c_p, Constraint_rod *c)
{
  int i;
  double fac, dist, vec[3], c_dist_2, c_dist;
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, c_p->r.type);

  /* also needed for coulomb */
  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = p1->r.p[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }
  vec[2] = 0.;

  /* charge stuff. This happens even if the particle does not feel the constraint. The electrostatic
     formulas for pbc highly dislike partial interactions anyways.
     THIS HAS TO BE DONE FIRST SINCE LJ CHANGES vec!!!
  */
  /*  fprintf(stderr, "%d: bj %f q %f l %f\n", this_node, coulomb.bjerrum, p1->r.q, c->lambda);*/
#ifdef ELECTROSTATICS
  if (coulomb.bjerrum > 0.0 && p1->r.q != 0.0 && c->lambda != 0.0) {
    fac = 2*coulomb.bjerrum*c->lambda*p1->r.q/c_dist_2;
    if (temperature > 0)
      fac *= temperature;
    p1->f[0]  += fac*vec[0];
    p1->f[1]  += fac*vec[1];
    c_p->f[0] -= fac*vec[0];
    c_p->f[1] -= fac*vec[1];
    /* fprintf(stderr, "%d: vec %f %f -> f %f %f\n", this_node, vec[0], vec[1],
       fac*vec[0], fac*vec[1]); */
  }
#endif
  if (ia_params->LJ_cut > 0. ) {
    /* put an infinite cylinder along z axis */
    c_dist = sqrt(c_dist_2);
    dist = c_dist - c->rad;

    if (dist > 0) {
      fac = dist / c_dist;
      for(i=0;i<2;i++) vec[i] *= fac;
      add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else {
      fprintf(stderr,"CONSTRAINT ROD: ERROR! part %d at (%.2e,%.2e,%.2e) within rod!\n",
	      p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
      errexit();
    }
  }
}

MDINLINE void add_constraints_forces(Particle *p1)
{
  int n=0;
  double dist, vec[3];
  IA_parameters *ia_params;
    
  for(n=0;n<n_constraints;n++) {
    ia_params=get_ia_param(p1->r.type, (&constraints[n].part_rep)->r.type);
    dist=0.;
	
    if(ia_params->LJ_cut > 0. ) {
      switch(constraints[n].type) {
      case CONSTRAINT_WAL: 
        calculate_wall_dist(p1, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
        if (dist > 0)
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist);
        else {
	  fprintf(stderr,"CONSTRAINT WALL : ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	  errexit();
        }
	break;
    
      case CONSTRAINT_SPH: 
        calculate_sphere_dist(p1, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
        if (dist > 0) {
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist);
        }
        else {
	  fprintf(stderr,"CONSTRAINT SPHERE: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	  errexit();
        }
	break;
    
      case CONSTRAINT_CYL: 
        calculate_cylinder_dist(p1, &constraints[n].part_rep, &constraints[n].c.cyl, &dist , vec); 
        if ( dist > 0 ) {
	  add_lj_pair_force(&constraints[n].part_rep, p1, ia_params, vec, dist);
        }
        else {
	  fprintf(stderr,"CONSTRAINT CYLINDER: ERROR! part %d at (%.2e,%.2e,%.2e) violated the constraint! dist= %f\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2],dist);
	  errexit();
        }
	break;

      case CONSTRAINT_MAZE: 
        calculate_maze_dist(p1, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
        if (dist > 0) {
	  add_lj_pair_force(p1, &constraints[n].part_rep, ia_params, vec, dist);
        }
        else {
	  fprintf(stderr,"CONSTRAINT MAZE: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	  errexit();
        }
	break;
      }
    }
    if (constraints[n].type == CONSTRAINT_ROD)
      add_rod_force(p1, &constraints[n].part_rep, &constraints[n].c.rod);
  }
}

MDINLINE double add_constraints_energy(Particle *p1, int n)
{
  double dist, vec[3];
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, (&constraints[n].part_rep)->r.type);
  dist=0.;
  if(ia_params->LJ_cut > 0. ) {
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      calculate_wall_dist(p1, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
      if (dist > 0)
	return lj_pair_energy(p1, &constraints[n].part_rep, ia_params, vec, dist);
      else {
	fprintf(stderr,"CONSTRAINT WALL : ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	errexit();
      }
      break;
    
    case CONSTRAINT_SPH: 
      calculate_sphere_dist(p1, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
      if (dist > 0) {
	return lj_pair_energy(p1, &constraints[n].part_rep, ia_params, vec, dist);
      }
      else {
	fprintf(stderr,"CONSTRAINT SPHERE: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	errexit();
      }
      break;
    
    case CONSTRAINT_CYL: 
      calculate_cylinder_dist(p1, &constraints[n].part_rep, &constraints[n].c.cyl, &dist , vec); 
      if ( dist > 0 ) {
	return lj_pair_energy(&constraints[n].part_rep, p1, ia_params, vec, dist);
      }
      else {
	fprintf(stderr,"CONSTRAINT CYLINDER: ERROR! part %d at (%.2e,%.2e,%.2e) violated the constraint! dist= %f\n",
		p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2],dist);
	errexit();
      }
      break;
	
    }
  }

  return 0.;
}

MDINLINE void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f[i] = 0;
}
#endif

#endif
