
// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#include "virtual_sites.h"
#include "pressure.h"

#ifdef VIRTUAL_SITES

// The following four functions are independent of the specif
// rules used to place virtual particles

void update_mol_vel_pos()
{
   //replace this by a better implementation later!
   update_mol_vel();
   update_mol_pos();
}

void update_mol_vel()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       update_mol_vel_particle(&p[i]);
    }
  }
}

void update_mol_pos()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       update_mol_pos_particle(&p[i]);
    }
    //only for real particles
  }
}

int update_mol_pos_cfg(){
  Particle *p;
  int i,j;
  double r_com[3];
  for(i=0; i<n_total_particles; i++) {
     p=&partCfg[i];
     update_mol_pos_particle(p);
  }
  return 1;
}

// Now, load the rules deciding, how to place virtual particles
// and transfer back forces to the real particles
#ifdef VIRTUAL_SITES_COM
 #include "virtual_sites_com.c"
#endif
#ifdef VIRTUAL_SITES_RELATIVE
 #include "virtual_sites_relative.c"
#endif

#endif 

