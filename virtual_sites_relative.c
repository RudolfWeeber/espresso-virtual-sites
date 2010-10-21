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
 return local_particles[p->p.mol_id];
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
  fprintf(stderr,"virtual_sites_relative.c - update_mol_pos_particle(): No real particle associated with virtual site.\n");
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
 
 // Calculate the new position of the virtual sites from
 // position of real particle + director 
 int i;
 for (i=0;i<3;i++)
  p->r.p[i] =p_real->r.p[i] +director[i];
}

// Update the vel of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_vel_particle(Particle *p)
{}

// Distribute forces that have accumulated on virtual particles to the 
// associated real particles
void distribute_mol_force()
{}

#endif

