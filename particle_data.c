/** \file particle_data.c
    This file contains everything related to particle storage. If you want to add a new
    property to the particles, it is probably a good idea to modify \ref part to give
    scripts access to that property. You always have to modify two positions: first the
    print section, where you should add your new data at the end, and second the read
    section where you have to find a nice and short name for your property to appear in
    the Tcl code. Then you just parse your part out of argc and argv.

    The corresponding header file is \ref particle_data.h "particle_data.h".
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "particle_data.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include "interaction_data.h" 
#include "debug.h"
#include "utils.h"
#include "cells.h"

/* cwz-build-command: make all
 */

/************************************************
 * defines
 ************************************************/

/** granularity of the particle buffers in particles */
#define PART_INCREMENT 32

/************************************************
 * variables
 ************************************************/

int max_seen_particle = -1;
int max_particle_node = 0;
int *particle_node = NULL;
Particle **local_particles = NULL;

/************************************************
 * functions
 ************************************************/

void map_particle_node(int part, int node)
{
  if (part > max_seen_particle) {
    int old_max = max_seen_particle, i;

    max_seen_particle = part;
    mpi_bcast_parameter(FIELD_MAXPART);

    if (part >= max_particle_node) {
      /* round up part + 1 in granularity PART_INCREMENT */
      max_particle_node = PART_INCREMENT*((part + PART_INCREMENT)/PART_INCREMENT);
      particle_node = (int *)realloc(particle_node, sizeof(int)*max_particle_node);
    }
    for (i = old_max + 1; i < max_seen_particle; i++)
      particle_node[i] = -1;
  }
  PART_TRACE(fprintf(stderr, "mapping %d -> %d (%d)\n", part, node, max_seen_particle));
  particle_node[part] = node;
}

void build_particle_node()
{
  if (max_seen_particle == -1)
    return;
  if (particle_node)
    free(particle_node);

  /* round up max_seen_particle + 1 in granularity PART_INCREMENT */
  max_particle_node = PART_INCREMENT*((max_seen_particle + PART_INCREMENT)/PART_INCREMENT);
  particle_node = malloc(max_particle_node*sizeof(int));
  mpi_who_has();
}

void init_particleList(ParticleList *pList)
{
  pList->n    = 0;
  pList->max  = 0;
  pList->part = NULL;
}

void realloc_particles(ParticleList *l, int size)
{
  int old_max = l->max, i;
  if (size < l->max) {
    /* shrink not as fast, just lose half, rounded up */
    l->max = PART_INCREMENT*(((l->n + size)/2 +
			      PART_INCREMENT - 1)/PART_INCREMENT);
  }
  else
    /* round up */
    l->max = PART_INCREMENT*((size + PART_INCREMENT - 1)/PART_INCREMENT);
  if (l->max != old_max)
    l->part = (Particle *) realloc(l->part, sizeof(Particle)*l->max);
  for (i = old_max; i < l->max; i++)
    l->part[i].r.identity = -1;
}

void init_redParticleList(RedParticleList *pList)
{
  pList->n    = 0;
  pList->max  = 0;
  pList->part = NULL;
}

void realloc_redParticles(RedParticleList *pList, int size)
{
  int old_max = pList->max, i;
  if (size < pList->max) {
    /* shrink not as fast, just lose half, rounded up */
    pList->max = PART_INCREMENT*(((pList->n + size)/2 +
			      PART_INCREMENT - 1)/PART_INCREMENT);
  }
  else
    /* round up */
    pList->max = PART_INCREMENT*((size + PART_INCREMENT - 1)/PART_INCREMENT);
  if (pList->max != old_max)
    pList->part = (ReducedParticle *) realloc(pList->part, sizeof(ReducedParticle)*pList->max);
  for (i = old_max; i < pList->max; i++)
    pList->part[i].identity = -1;
}

int try_delete_bond(Particle *part, int *bond)
{
  IntList *bl = &part->bl;
  int i, j, type, partners;
  for (i = 0; i < bl->n;) {
    type = bond[i];
    partners = bonded_ia_params[type].num;
    if (type != bond[0])
      i += 1 + partners;
    else {
      for(j = 1; j <= partners; j++)
	if (bond[j] != bl->e[i + j])
	  break;
      if (j > partners) {
	bl->n -= 1 + partners;
	memcpy(bl->e + i, bl->e + i + 1 + partners,
	       bl->n - i);
	realloc_intlist(bl, bl->n);
	return 1;
      }
      i += 1 + partners;
    }
  }
  return 0;
}

Particle *got_particle(ParticleList *l, int id)
{
  int i;

  for (i = 0; i < l->n; i++)
    if (l->part[i].r.identity == id)
      break;
  if (i == l->n)
    return NULL;
  return &(l->part[i]);
}

Particle *alloc_particle(ParticleList *l)
{
  int index;

  index = l->n++;
  realloc_particles(l, l->n);

  return &l->part[index];
}

Particle *append_particle(ParticleList *l, Particle *part)
{
  Particle *p = alloc_particle(l);
  memcpy(p, part, sizeof(Particle));
  return p;
}

Particle *move_particle(ParticleList *dl, ParticleList *sl, int i)
{
  Particle *dst = alloc_particle(dl);
  Particle *src = &sl->part[i];
  Particle *end = &sl->part[sl->n - 1];
  memcpy(dst, src, sizeof(Particle));
  if (--sl->n > 0)
    memcpy(src, end, sizeof(Particle));
  realloc_particles(sl, sl->n);
  return dst;
}

void fold_particle(double pos[3],int image_box[3])
{
  int i;
  int tmp;
  for(i=0;i<3;i++) {
    if (periodic[i]) {
      image_box[i] += (tmp = (int)floor(pos[i]/box_l[i]));
      pos[i]       = pos[i] - tmp*box_l[i];    
      if(pos[i] < 0 || pos[i] > box_l[i]) {
	fprintf(stderr,"%d: fold_particle: Particle out of"
		" range image_box[%d] = %d, exiting\n",
		this_node,i,image_box[i]);
	errexit();
      }
    }
  }
}

void unfold_particle(double pos[3],int image_box[3])
{
  int i;
  for(i=0;i<3;i++) {
    pos[i]       = pos[i] + image_box[i]*box_l[i];    
    image_box[i] = 0;
  }
}

void particle_finalize_data()
{
  /* happens if 0 particles... */
  if (!node_grid_is_set())
    setup_node_grid();

  mpi_bcast_parameter(FIELD_MAXPART);

  /* invalidate particle->node data */
  if (particle_node) {
    free(particle_node);
    particle_node = NULL;
  }
}

void local_particles_init()
{
  Particle *p;
  int np, m, n, o, i;
  if (max_seen_particle >= 0)
    local_particles = (Particle **)malloc((max_seen_particle + 1)*sizeof(Particle *));
  for(i = 0; i <= max_seen_particle; i++) local_particles[i] = NULL;
  INNER_CELLS_LOOP(m, n, o) {
    p = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) local_particles[p[i].r.identity] = & p[i];
  }
}

void local_particles_exit()
{
  if (local_particles) {
    free(local_particles);
    local_particles = NULL;
  }
}

/** append particle data in ASCII form to the Tcl result.
    @param part_num the particle which data is appended
    @param interp   the Tcl interpreter to which result to add to */
int printParticleToResult(Tcl_Interp *interp, int part_num)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle part;
  IntList *bl = &(part.bl);
  int node;

  if (part_num < 0 || part_num > max_seen_particle)
    return (TCL_ERROR);

  if (!particle_node)
    build_particle_node();

  node = particle_node[part_num];
  if (node == -1)
    return (TCL_ERROR);

  mpi_recv_part(node, part_num, &part);
  sprintf(buffer, "%d", part.r.identity);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  Tcl_PrintDouble(interp, part.r.p[0], buffer);
  Tcl_AppendResult(interp, " pos ", buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.p[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.p[2], buffer);
  Tcl_AppendResult(interp, buffer, " type ", (char *)NULL);
  sprintf(buffer, "%d", part.r.type);
  Tcl_AppendResult(interp, buffer, " q ", (char *)NULL);
  Tcl_PrintDouble(interp, part.r.q, buffer);
  Tcl_AppendResult(interp, buffer, " v ", (char *)NULL);
  Tcl_PrintDouble(interp, part.v[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.v[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.v[2], buffer);
  Tcl_AppendResult(interp, buffer, " f ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, part.f[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  /* print bonding structure */
  if(bl->n > 0) {
    int i=0,j,size;
    Tcl_AppendResult(interp, " bonds { ", (char *)NULL);
    while(i<bl->n) {
      size = bonded_ia_params[bl->e[i]].num;
      sprintf(buffer, "{%d ", bl->e[i]); i++;
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      for(j=0;j<size-1;j++) {
	sprintf(buffer, "%d ", bl->e[i]); i++;
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      sprintf(buffer, "%d} ", bl->e[i]); i++;
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
  realloc_intlist(bl, 0);
  return (TCL_OK);
}

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;
  int node, j;

  if (!node_grid_is_set())
    setup_node_grid();

  if (!particle_node)
    build_particle_node();

  /* if no further arguments are given, print out all stored particles */
  if (argc == 1) {
    int i = 0, start = 1;
    PART_TRACE(fprintf(stderr, "max_seen %d\n", max_seen_particle));
    for (i = 0; i <= max_seen_particle ; i++) {
      PART_TRACE(fprintf(stderr, "particle %d\n", i));
      if (particle_node[i] != -1) {
	if (start) {
	  Tcl_AppendResult(interp, "{", (char *)NULL);
	  start = 0;
	}
	else
	  Tcl_AppendResult(interp, " {", (char *)NULL);
	printParticleToResult(interp, i);
	Tcl_AppendResult(interp, "}", (char *)NULL);
      }
    }
    return (TCL_OK);
  }

  if (argc == 2 && !strncmp(argv[1], "number", strlen(argv[1]))) {
    char buffer[TCL_INTEGER_SPACE];
    int i = 0, c = 0;
    for (i = 0; i <= max_seen_particle; i++)
      if (particle_node[i] != -1)
	c++;
    sprintf(buffer, "%d", c);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }

  part_num = atol(argv[1]);

  /* print out particle information */
  if (argc == 2) {
    if (printParticleToResult(interp, part_num) == TCL_ERROR)
      Tcl_AppendResult(interp, "na", (char *)NULL);
    return (TCL_OK);
  }

  /* print out partial particle information */
  if (!strncmp(argv[2], "print",  strlen(argv[2]))) {
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Particle part;
    IntList *bl = &(part.bl);

    if (part_num < 0 || part_num > max_seen_particle) {
      Tcl_AppendResult(interp, "na", (char *)NULL);
      return (TCL_OK);
    }

    if (!particle_node)
      build_particle_node();

    node = particle_node[part_num];
    if (node == -1) {
      Tcl_AppendResult(interp, "na", (char *)NULL);
      return (TCL_OK);
    }

    mpi_recv_part(node, part_num, &part);

    argc -= 3;
    argv += 3;
    while (argc > 0) {
      if (!strncmp(argv[0], "identity", strlen(argv[0]))) {
	sprintf(buffer, "%d", part.r.identity);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "pos", strlen(argv[0]))) {
	Tcl_PrintDouble(interp, part.r.p[0], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.r.p[1], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.r.p[2], buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "type", strlen(argv[0]))) {
	sprintf(buffer, "%d", part.r.type);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "q", strlen(argv[0]))) {
	Tcl_PrintDouble(interp, part.r.q, buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "v", strlen(argv[0]))) {
	Tcl_PrintDouble(interp, part.v[0], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.v[1], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.v[2], buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "f", strlen(argv[0]))) {
	Tcl_PrintDouble(interp, part.f[0], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.f[1], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	Tcl_PrintDouble(interp, part.f[2], buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      else if (!strncmp(argv[0], "bonds", strlen(argv[0]))) {
	int i = 0, j, size;
	Tcl_AppendResult(interp, "{", (char *)NULL);
	while(i < bl->n) {
	  size = bonded_ia_params[bl->e[i]].num;
	  sprintf(buffer, "{%d ", bl->e[i]);
	  i++;
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	  for(j = 0; j < size - 1; j++) {
	    sprintf(buffer, "%d ", bl->e[i]);
	    i++;
	    Tcl_AppendResult(interp, buffer, (char *)NULL);
	  }
	  sprintf(buffer, "%d} ", bl->e[i]);
	  i++;
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	}
	Tcl_AppendResult(interp, "}", (char *)NULL);
	realloc_intlist(bl, 0);
      }
      else {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "unknown particle data \"", argv[0], "\" requested", (char *)NULL);
      }
      if (argc > 1)
	Tcl_AppendResult(interp, " ", (char *)NULL);

      argc--;
      argv++;
    }
  }

  /* set particle data */
  if (part_num < 0) {
    Tcl_AppendResult(interp, "particle identities must be positive", (char *)NULL);
    return (TCL_ERROR);
  }
  node = (part_num <= max_seen_particle) ? particle_node[part_num] : -1;

  argc -= 2;
  argv += 2;
  while (argc > 0) {
    if (!strncmp(argv[0], "pos", strlen(argv[0]))) {
      double pos[3];
      if (argc < 4) {
	Tcl_AppendResult(interp, "pos requires 3 arguments", (char *) NULL);
	return (TCL_ERROR);
      }
      /* set position */
      for (j = 0; j < 3; j++) {
	if (Tcl_GetDouble(interp, argv[1 + j], &pos[j]) == TCL_ERROR)
	  return (TCL_ERROR);
      }
    
      if (node == -1) {
	/* spatial position */
	node = find_node(pos);
	map_particle_node(part_num, node);
      }
      mpi_place_particle(node, part_num, pos);

      argc -= 4;
      argv += 4;
    }
    else {
      if (node == -1) {
	Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
	return (TCL_ERROR);
      }
      
      if (!strncmp(argv[0], "q", strlen(argv[0]))) {
	double q;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "q requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set charge */
	if (Tcl_GetDouble(interp, argv[1], &q) == TCL_ERROR)
	  return (TCL_ERROR);
	
	mpi_send_q(node, part_num, q);

	argc -= 2;
	argv += 2;
      }
      else if (!strncmp(argv[0], "v", strlen(argv[0]))) {
	double v[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "v requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &v[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &v[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &v[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_v(node, part_num, v);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "f", strlen(argv[0]))) {
	double f[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "f requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &f[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &f[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &f[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_f(node, part_num, f);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "type", strlen(argv[0]))) {
	int type;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "type requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set type */
	if (Tcl_GetInt(interp, argv[1], &type) == TCL_ERROR)
	  return (TCL_ERROR);

	if (type < 0) {
	  Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
	  return (TCL_ERROR);	  
	} 

	/* make sure type exists */
	make_particle_type_exist(type);

	mpi_send_type(node, part_num, type);

	argc -= 2;
	argv += 2;
      }
      else if (!strncmp(argv[0], "bond", strlen(argv[0]))) {
	int delete = 0;
	int type_num;
	int n_partners;
	/* Bond type number and the bond partner atoms are stored in this field. */
	int *bond;
	/* check particle position */
	if (node == -1) {
	  Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
	  return (TCL_ERROR);
	}
	/* check number of arguments */
	if (argc < 3) {
	  Tcl_AppendResult(interp, "bond requires at least 2 arguments: "
			   "[delete] <type_num> <partner>", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* parse away delete eventually */
	if (!strncmp(argv[1], "delete", strlen(argv[1]))) {
	  delete = 1;
	  argc--;
	  argv++;
	}
	/* check type_num */
	if (Tcl_GetInt(interp, argv[1], &type_num) == TCL_ERROR)
	  return (TCL_ERROR);
	if(type_num < 0 || type_num >= n_bonded_ia) {
	  Tcl_AppendResult(interp, "invalid bonded interaction type_num"
			   "(Set bonded interaction parameters first)", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* check partners */ 
	n_partners = bonded_ia_params[type_num].num;
	if(argc < 2+n_partners) {
	  char buffer[256 + 2*TCL_INTEGER_SPACE];
	  sprintf(buffer, "bond type %d requires %d arguments.",
		  type_num, n_partners+1);
	  Tcl_AppendResult(interp, buffer, (char *) NULL);
	  return (TCL_ERROR);
	}
	bond = (int *)malloc( (n_partners+1)*sizeof(int) );
	bond[0] = type_num;
	j=1;
	while(j <= n_partners) {
	  if (Tcl_GetInt(interp, argv[j+1], &(bond[j])) == TCL_ERROR) {
	    free(bond);
	    return (TCL_ERROR);
	  }
	  if(bond[j] > max_seen_particle || particle_node[bond[j]] == -1) {
	    Tcl_AppendResult(interp, "partner atom %d (identity %d) not known"
			     ,j+1,bond[j],
			     "(Set all partner atoms first)", (char *) NULL);
	    free(bond);
	    return (TCL_ERROR);
	  }
	  j++;
	}
	/* set/delete bond */ 
	if (!mpi_send_bond(node, part_num, bond, delete)) {
	  Tcl_AppendResult(interp, "bond to delete did not exist", (char *)NULL);
	  free(bond);
	  return (TCL_ERROR);
	}
	free(bond);
	argc -= (2 + n_partners);
	argv += (2 + n_partners);
      }
      else {
	Tcl_AppendResult(interp, "unknown particle parameter \"", argv[0],"\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }

  return (TCL_OK);
}
