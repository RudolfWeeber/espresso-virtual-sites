// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#ifdef VIRTUAL_SITES_RELATIVE

#include "rotation.h"

// This is the "relative" implementation for virtual sites.
// Virtual particles are placed relative to the position of a real particle

// Obtain the real particle from which a virtual particle derives it's position
// Note: for now, we use the mol_di property of Particle
Particle* vs_relative_get_real_particle(Particle* p)
{
 return local_particles[p->p.vs_relative_to_particle_id];
}




// Update the pos of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_pos_particle(Particle *p)
{
 // First obtain the real particle responsible for this virtual particle:
 // Find the 1st real particle in the topology for the virtual particle's mol_id
 Particle *p_real = vs_relative_get_real_particle(p);
 // Check, if a real particle was found
 if (!p_real)
 {
  printf(stderr,"virtual_sites_relative.c - update_mol_pos_particle(): No real particle associated with virtual site.\n");
  return;
 }
 
 // Calculate the quaternion defining the orientation of the vecotr connectinhg
 // the virtual site and the real particle
 // This is obtained, by multiplying the quaternion representing the director
 // of the real particle with the quaternion of the virtual particle, which 
 // specifies the relative orientation.
 double q[4];
 multiply_quaternions(p_real->r.quat,p->r.quat,q);
 // Calculate the director resulting from the quaternions
 double director[3];
 convert_quat_to_quatu(q,director);
 // normalize
 double l =sqrt(sqrlen(director));
 // Division comes in the loop below
 
 // Calculate the new position of the virtual sites from
 // position of real particle + director 
 int i;
 double new_pos[3];
 double tmp;
 for (i=0;i<3;i++)
 {
  new_pos[i] =p_real->r.p[i] +director[i]/l*p->p.vs_relative_distance;
  // Handle the case that one of the particles had gone over the periodic
  // boundary and its coordinate has been folded
  #ifdef PARTIAL_PERIODIC
   if (PERIODIC(i)) 
  #endif
  {
    tmp =p->r.p[i] -new_pos[i];
    //printf("%f\n",tmp);
    if (tmp > box_l[i]/2.) {
     //printf("greater than box_l/2 %f\n",tmp);
     p->r.p[i] =new_pos[i] + box_l[i];
    }
    else if (tmp < -box_l[i]/2.) {
     //printf("smaller than box_l/2 %f\n",tmp);
     p->r.p[i] =new_pos[i] - box_l[i];
    }
    else p->r.p[i] =new_pos[i];
   }

 }
}

// Update the vel of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_vel_particle(Particle *p)
{
 // NOT TESTED!
 
 // First obtain the real particle responsible for this virtual particle:
 Particle *p_real = vs_relative_get_real_particle(p);
 // Check, if a real particle was found
 if (!p_real)
 {
  printf(stderr,"virtual_sites_relative.c - update_mol_pos_particle(): No real particle associated with virtual site.\n");
  return;
 }
 
 // Calculate the quaternion defining the orientation of the vecotr connectinhg
 // the virtual site and the real particle
 // This is obtained, by multiplying the quaternion representing the director
 // of the real particle with the quaternion of the virtual particle, which 
 // specifies the relative orientation.
 double q[4];
 multiply_quaternions(p_real->r.quat,p->r.quat,q);
 // Calculate the director resulting from the quaternions
 double director[3];
 convert_quat_to_quatu(q,director);
 // normalize
 double l =sqrt(sqrlen(director));
 // Division comes in the loop below

 // Obtain velocity from v=v_real particle + omega_real_particle \times director
 vector_product(p_real->m.omega,director,p->m.v);

 int i;
 // Add prefactors and add velocity of real particle
 for (i=0;i<3;i++)
 {
  // Scale the velocity by the distance of virtual particle from the real particle
  // Also, espresso stores not velocity but velocity * time_step
  p->m.v[i] *= time_step * p->p.vs_relative_distance/l;
  // Add velocity of real particle
  p->m.v[i] += p_real->m.v[i];
 }
}


// Distribute forces that have accumulated on virtual particles to the 
// associated real particles
void distribute_mol_force()
{
  // Iterate over all the particles in the local cells
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      // We only care about virtual particles
      if (ifParticleIsVirtual(&p[i])) {

       // First obtain the real particle responsible for this virtual particle:
       Particle *p_real = vs_relative_get_real_particle(&p[i]);

       // Get distance vector pointing from real to virtual particle, respecting periodic boundary i
       // conditions
       double d[3];
       get_mi_vector(d,p[i].r.p,p_real->r.p);

       // The rules for transfering forces are:
       // F_realParticle +=F_virtualParticle
       // T_realParticle +=f_realParticle \times (r_virtualParticle-r_realParticle)
       
       // Calculate torque to be added on real particle
       double tmp[3];
       vector_product(d,p[i].f.f,tmp);

       // Add forces and torques
       int j;
       for (j=0;j<3;j++) {
         p_real->f.torque[j]+=tmp[j];
	 p_real->f.f[j]+=p[i].f.f[j];
	 // Clear forces on virtual particle
	 p[i].f.f[j]=0;

       }
      }
    }
  }
}



// Setup the virtual_sites_relative properties of a particle so that the given virtaul particle will follow the given real particle
int vs_relate_to(int part_num, int relate_to)
{
    // Get the data for the particle we act on and the one we wnat to relate
    // it to.
    Particle  p_current,p_relate_to;
    if ((get_particle_data(relate_to,&p_relate_to)!=TCL_OK) || 
        (get_particle_data(part_num,&p_current)!=TCL_OK)) {
         printf("Could not retrieve particle data for the given id.\n");
          return TCL_ERROR;
    }
    
    // get teh distance between the particles
    double d[3];
    get_mi_vector(d, p_current.r.p,p_relate_to.r.p);
    
    // Set the particle id of the particle we want to relate to and the distnace
    if (set_particle_vs_relative(part_num, relate_to, sqrt(sqrlen(d))) == TCL_ERROR) {
      printf("setting the vs_relative attributes failed.\n");
      return TCL_ERROR;
    }

    // Now, calculate the quaternions which specify the angle between 
    // the director of the particel we relate to and the vector
    // (paritlce_we_relate_to - this_particle)
    double quat[4];
    // The vs_relative implemnation later obtains the direcotr by multiplying
    // the quaternions representing the orientation of the real particle
    // with those in the virtual particle. The re quulting quaternion is then
    // converted to a director.
    // Whe have quat_(real particle) *quat_(virtual particle) 
    // = quat_(obtained from desired director)
    // Resolving this for the quat_(virtaul particle)

    //Normalize desired director
    double l=sqrt(sqrlen(d));
    int i;
    for (i=0;i<3;i++)
     d[i]/=l;

    // Obtain quaternions from desired director
    double quat_director[4];
    convert_quatu_to_quat(d, quat_director);

    // Define quat as described above:
    double x=0;
    for (i=0;i<4;i++)
     x+=p_relate_to.r.quat[i]*p_relate_to.r.quat[i];

    quat[0]=0;
    for (i=0;i<4;i++)
     quat[0] +=p_relate_to.r.quat[i]*quat_director[i];
    
    quat[1] =-quat_director[0] *p_relate_to.r.quat[1] 
       +quat_director[1] *p_relate_to.r.quat[0]
       +quat_director[2] *p_relate_to.r.quat[3]
       -quat_director[3] *p_relate_to.r.quat[2];
    quat[2] =p_relate_to.r.quat[1] *quat_director[3] 
      + p_relate_to.r.quat[0] *quat_director[2] 
      - p_relate_to.r.quat[3] *quat_director[1] 
      - p_relate_to.r.quat[2] * quat_director[0];
    quat[3] =quat_director[3] *p_relate_to.r.quat[0]
      - p_relate_to.r.quat[3] *quat_director[0] 
      + p_relate_to.r.quat[2] * quat_director[1] 
      - p_relate_to.r.quat[1] *quat_director[2];
    for (i=0;i<4;i++)
     quat[i]/=x;
   
   
   // Verify result
   double qtemp[4];
   multiply_quaternions(p_relate_to.r.quat,quat,qtemp);
   for (i=0;i<4;i++)
    if (fabs(qtemp[i]-quat_director[i])>1E-9)
     printf("Component %d: %f instead of %f\n",i,qtemp[i],quat_director[i]);


   // Save the quaternions in the particle
   if (set_particle_quat(part_num, quat) == TCL_ERROR) {
     printf("set particle position first\n");

     return TCL_ERROR;
   }
   
   return TCL_OK;
}

#endif

