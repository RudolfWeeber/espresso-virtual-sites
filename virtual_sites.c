// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#include "virtual_sites.h"
#include "pressure.h"

/** \file virtual_sites.h
 *  This file contains routine to handle virtual sites
 *  Virtual sites are like particles, but they will be not integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move no-virtual particles
 *  - update virtual sites
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:junghans@mpip-mainz.mpg.de">Ben</a>
 */

/** \name Private Functions */
void calc_mol_vel_pos(int mol_id,double r_com[3],double v_com[3]);
void calc_mol_vel(int mol_id,double v_com[3]);
void calc_mol_pos(int mol_id,double r_com[3]);
void put_mol_force_on_parts(Particle *p_com);

Particle *get_mol_com_particle(Particle *calling_p);
int get_mol_com_id(Particle *calling_p);


#ifdef VIRTUAL_SITES
void update_mol_vel_pos()
{
  Particle *p;
  int i, np, c, mol_id;
  Cell *cell;
  int j;
  double r_com[3],v_com[3];

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         mol_id = p[i].p.mol_id;
         calc_mol_vel_pos(mol_id,r_com,v_com);
         for (j=0;j<3;j++){
            p[i].r.p[j] = r_com[j];
            p[i].m.v[j] = v_com[j];
         }
      }
    }
  }
}

void update_mol_vel()
{
  Particle *p;
  int i, np, c, mol_id;
  Cell *cell;
  int j;
  double v_com[3];

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         mol_id = p[i].p.mol_id;
         calc_mol_vel(mol_id,v_com);
         for (j=0;j<3;j++){
            p[i].m.v[j] = v_com[j];
         }
      }
    }
  }
}

void update_mol_pos()
{
  Particle *p;
  int i, np, c, mol_id;
  Cell *cell;
  int j;
  double r_com[3];

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         mol_id = p[i].p.mol_id;
         calc_mol_pos(mol_id,r_com);
         for (j=0;j<3;j++){
            p[i].r.p[j] = r_com[j];
         }
      }
    }
  }
}

void distribute_mol_force()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         if (sqrlen(p[i].f.f)!=0){
            put_mol_force_on_parts(&p[i]);
         }
      }
    }
  }
}

void calc_mol_vel_pos(int mol_id,double r_com[3],double v_com[3]){
   int i,j;
   double M=0;
   Particle *p;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
   for (i=0;i<3;i++){
      r_com[i]=0.0;
      v_com[i]=0.0;
   }
   for (i=0;i<topology[mol_id].part.n;i++){
       p=local_particles[topology[mol_id].part.e[i]];
       if (!ifParticleIsVirtual(p)) {
          for (j=0;j<3;j++){
              unfold_position(p->r.p,p->l.i);
              r_com[j] += PMASS(*p)*p->r.p[j];
              fold_position(p->r.p,p->l.i);
              v_com[j] += PMASS(*p)*p->m.v[j];
          }
          M+=PMASS(*p);
#ifdef VIRTUAL_SITES_DEBUG
          count++;
#endif
       }
   }
   for (j=0;j<3;j++){
      r_com[j] /= M;
      v_com[j] /= M;
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
      fprintf(stderr,"There is more than one COM in calc_mol_vel_pos! mol_id=%i\n",mol_id);
      exit(182);
   }
#endif
}

void calc_mol_vel(int mol_id,double v_com[3]){
   int i,j;
   double M=0;
   Particle *p;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
   for (i=0;i<3;i++){
      v_com[i]=0.0;
   }
   for (i=0;i<topology[mol_id].part.n;i++){
       p=local_particles[topology[mol_id].part.e[i]];
       if (!ifParticleIsVirtual(p)) {
          for (j=0;j<3;j++){
              v_com[j] += PMASS(*p)*p->m.v[j];
          }
          M+=PMASS(*p);
#ifdef VIRTUAL_SITES_DEBUG
          count++;
#endif
       }
   }
   for (j=0;j<3;j++){
      v_com[j] /= M;
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
      fprintf(stderr,"There is more than one COM in calc_mol_vel! mol_id=%i\n",mol_id);
      exit(182);
   }
#endif
}

void calc_mol_pos(int mol_id,double r_com[3]){
   int i,j;
   double M=0;
   Particle *p;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
   for (i=0;i<3;i++){
      r_com[i]=0.0;
   }
   for (i=0;i<topology[mol_id].part.n;i++){
       p=local_particles[topology[mol_id].part.e[i]];
       if (!ifParticleIsVirtual(p)) {
          for (j=0;j<3;j++){
              unfold_position(p->r.p,p->l.i);
              r_com[j] += PMASS(*p)*p->r.p[j];
              fold_position(p->r.p,p->l.i);
          }
          M+=PMASS(*p);
#ifdef VIRTUAL_SITES_DEBUG
          count++;
#endif
       }
   }
   for (j=0;j<3;j++){
      r_com[j] /= M;
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
      fprintf(stderr,"There is more than one COM in calc_mol_pos! mol_id=%i\n",mol_id);
      exit(182);
   }
#endif
}

void put_mol_force_on_parts(Particle *p_com){
   int i,j,mol_id;
   Particle *p;
   double force[3],fac;
 #ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
  mol_id=p_com->p.mol_id;
   fac=topology[mol_id].part.n-1;
   for (i=0;i<3;i++){
      force[i]=fac*p_com->f.f[i];
      p_com->f.f[i]=0.0;
   }
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];
      if (!ifParticleIsVirtual(p)) {
         for (j=0;j<3;j++){
            p->f.f[j]+=force[j];
         }
#ifdef VIRTUAL_SITES_DEBUG
         count++;
#endif
      }
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
      fprintf(stderr,"There is more than one COM input_mol_force_on_parts! mol_id=%i\n",mol_id);
      exit(182);
   }
#endif
}

Particle *get_mol_com_particle(Particle *calling_p){
   int mol_id;
   int i;
   Particle *p;
   mol_id=calling_p->p.mol_id;
   for (i=0;i<topology[mol_id].part.n;i++){
       p=local_particles[topology[mol_id].part.e[i]];
       if (ifParticleIsVirtual(p)) {
          return p;
       }
   }
#ifdef VIRTUAL_SITES_DEBUG
   fprintf(stderr,"No com found in get_mol_com ! pnr=%i\n",calling_p->p.identity);
   exit(182);
#endif
   return calling_p;
}

int get_mol_com_id(Particle *calling_p){
   Particle *p_com;
   int id;
   p_com=get_mol_com_particle(calling_p);
   id=p_com->p.identity;
   return id;
}

void get_mol_dist_vector(Particle *p1,Particle *p2,double dist[3]){
   Particle *p1_com,*p2_com;
   p1_com=get_mol_com_particle(p1);
   p2_com=get_mol_com_particle(p2);
   get_mi_vector(dist,p1_com->r.p, p2_com->r.p);
   //vecsub(p1_com->r.p,p2_com->r.p,dist);
}

double get_mol_dist(Particle *p1,Particle *p2){
   double dist[3],dist2;
   get_mol_dist_vector(p1,p2,dist);
   dist2=SQR(dist[0])+SQR(dist[1])+SQR(dist[2]);
   return sqrt(dist2);
}

/** \name Statistic Functions */
/** \name Private Functions */
double calc_pressure_mol(int type1,int type2);
double calc_energy_kinetic_mol(int type);
void calc_force_between_mol(int mol_id1,int mol_id2,double force[3]);

Particle *get_mol_com_particle_from_molid_cfg(int mol_id);
Particle *get_mol_com_particle_cfg(Particle *calling_p);
int get_mol_com_id_cfg(Particle *calling_p);
void get_mol_dist_vector_from_molid_cfg(int mol_id1,int mol_id2,double dist[3]);
void get_mol_dist_vector_cfg(Particle *p1,Particle *p2,double dist[3]);
double get_mol_dist_cfg(Particle *p1,Particle *p2);

int parse_and_print_pressure_mol(Tcl_Interp *interp,int argc, char **argv)
{
   char buffer[TCL_DOUBLE_SPACE];
   int type1, type2;
   double psum;
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_pressure_mol: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   if (argc < 2) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }

   if (!ARG0_IS_I(type1)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }
   if (!ARG1_IS_I(type2)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze pressure_mol <type1> <type2>", (char *)NULL);
      return (TCL_ERROR);
   }
   argc-=2; argv+=2;

   if (n_molecules==0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "No molecules defined !", (char *)NULL);
      return (TCL_ERROR);
   }
   psum=calc_pressure_mol(type1,type2);
   //sprintf(buffer,"%i",type1);
   //Tcl_AppendResult(interp,"{ analyze pressure_mol ",buffer," ",(char *)NULL);   
   //sprintf(buffer,"%i",type2);
   //Tcl_AppendResult(interp,buffer," ",(char *)NULL);   
   sprintf(buffer,"%e",psum);
   Tcl_AppendResult(interp, buffer,(char *)NULL);
   return TCL_OK;
}

int parse_and_print_energy_kinetic_mol(Tcl_Interp *interp,int argc, char **argv)
{
   char buffer[TCL_DOUBLE_SPACE];
   int type;
   double Ekin;
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_energy_kinetic_mol: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   if (argc < 1) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type>", (char *)NULL);
      return (TCL_ERROR);
   }

   if (!ARG0_IS_I(type)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type>", (char *)NULL);
      return (TCL_ERROR);
   }
   argc-=1; argv+=1;

   if (n_molecules==0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "No molecules defined !", (char *)NULL);
      return (TCL_ERROR);
   }
   Ekin=calc_energy_kinetic_mol(type);
   //sprintf(buffer,"%i",type);
   //Tcl_AppendResult(interp,"{ analyze pressure_mol ",buffer," ",(char *)NULL);   
   sprintf(buffer,"%e",Ekin);
   Tcl_AppendResult(interp, buffer,(char *)NULL);
   return TCL_OK;
}

double calc_pressure_mol(int type1,int type2){
   double force[3],com_dist[3],psum=0;
   int i,j,k,start;

   for (i=0;i<n_molecules;i++){
      if (topology[i].type == type1){
         if (type1==type2){
            start=i+1;
         } else {
            start=0;
         }
         for (j=start;j<n_molecules;j++){
            if (topology[j].type == type2){
               get_mol_dist_vector_from_molid_cfg(i,j,com_dist);
               calc_force_between_mol(i,j,force);
               for (k=0;k<3;k++){
                   psum+=force[k]*com_dist[k];
               }
            }
         }
      }
   }
   psum/=3.0;
   return psum;
}

void calc_force_between_mol(int mol_id1,int mol_id2,double force[3]){
   int i,j;
   Particle *p1,*p2;
   double vec12[3],dist2,dist;
   #ifdef ELECTROSTATICS
   #ifndef INTER_RF
   fprintf(stderr,"parse_and_print_pressure_mol is only possible with INTER_RF ");
   exit(182);
   #endif
   #endif

   force[0]=force[1]=force[2]=0.0;
   for (i=0;i<topology[mol_id1].part.n;i++){
      p1=&partCfg[topology[mol_id1].part.e[i]];
      for (j=0;j<topology[mol_id2].part.n;j++){
         p2=&partCfg[topology[mol_id2].part.e[j]];
         get_mi_vector(vec12,p1->r.p, p2->r.p);
         dist2=sqrlen(vec12);
         dist=sqrt(dist2);
         calc_non_bonded_pair_force_simple(p1,p2,vec12,dist,dist2,force);
      }
   }
}

double calc_energy_kinetic_mol(int type){
   double E_kin=0;
   int i;
   Particle *p_com;
   for (i=0;i<n_molecules;i++){
      if (topology[i].type == type){
         p_com=get_mol_com_particle_from_molid_cfg(i);
#ifdef VIRTUAL_SITES_DEBUG
         if (!ifParticleIsVirtual(p_com)){
            fprintf(stderr,"Could not fetch com in get_mol_dist_cfg! From %i\n",p_com->p.identity);
            exit(182);
         }
#endif
         E_kin+=PMASS(*p_com)*sqrlen(p_com->m.v);
      }
   }
   E_kin*=0.5/time_step/time_step;
   return E_kin;
}

Particle *get_mol_com_particle_from_molid_cfg(int mol_id){
   Particle *calling_p,*p_com;
   calling_p=&partCfg[topology[mol_id].part.e[0]];
   p_com=get_mol_com_particle_cfg(calling_p);
   return p_com;
}

Particle *get_mol_com_particle_cfg(Particle *calling_p){
   int mol_id;
   int i;
   Particle *p;
   mol_id=calling_p->p.mol_id;
   for (i=0;i<topology[mol_id].part.n;i++){
       p=&partCfg[topology[mol_id].part.e[i]];
       if (ifParticleIsVirtual(p)){
          return p;
       }
   }
#ifdef VIRTUAL_SITES_DEBUG
   fprintf(stderr,"No com found in get_mol_com ! pnr=%i\n",calling_p->p.identity);
   exit(182);
#endif
   return calling_p;
}

int get_mol_com_id_cfg(Particle *calling_p){
   Particle *p_com;
   int id;
   p_com=get_mol_com_particle_cfg(calling_p);
   id=p_com->p.identity;
   return id;
}

void get_mol_dist_vector_from_molid_cfg(int mol_id1,int mol_id2,double dist[3]){
   Particle *p1,*p2;
   p1=&partCfg[topology[mol_id1].part.e[0]];
   p2=&partCfg[topology[mol_id2].part.e[0]];
   get_mol_dist_vector_cfg(p1,p2,dist);
}

void get_mol_dist_vector_cfg(Particle *p1,Particle *p2,double dist[3]){
   Particle *p1_com,*p2_com;
   p1_com=get_mol_com_particle_cfg(p1);
#ifdef VIRTUAL_SITES_DEBUG
   if (!ifParticleIsVirtual(p1_com)){
      fprintf(stderr,"Could not fetch com in get_mol_dist_cfg! From %i to %i\n",p1->p.identity,p1_com->p.identity);
      exit(182);
   }
#endif
   p2_com=get_mol_com_particle_cfg(p2);
#ifdef VIRTUAL_SITES_DEBUG
   if (!ifParticleIsVirtual(p2_com)){
      fprintf(stderr,"Could not fetch com in get_mol_dist_cfg! From %i to %i\n",p2->p.identity,p2_com->p.identity);
      exit(182);
   }
#endif
   get_mi_vector(dist,p1_com->r.p, p2_com->r.p);
}

double get_mol_dist_cfg(Particle *p1,Particle *p2){
   double dist[3],dist2;
   get_mol_dist_vector(p1,p2,dist);
   dist2=SQR(dist[0])+SQR(dist[1])+SQR(dist[2]);
   return sqrt(dist2);
}

#endif
