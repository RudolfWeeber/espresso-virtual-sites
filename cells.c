// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file cells.c
 *
 *  This file contains everything related to the link cell
 *  algorithm. 
 *
 *  For more information on cells,
 *  see \ref cells.h "cells.h"
 *   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cells.h"
#include "config.h"
#include "debug.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "communication.h"
#include "utils.h"
#include "verlet.h"
#include "ghosts.h"
#include "parser.h"
#include "domain_decomposition.h"
#include "nsquare.h"

/* Variables */

/** list of all cells. */
Cell *cells = NULL;
/** size of \ref cells */
int n_cells = 0;
/** list of pointers to all cells containing particles physically on the local node. */
CellPList local_cells = { NULL, 0, 0 };
/** list of pointers to all cells containing ghosts. */
CellPList ghost_cells = { NULL, 0, 0 };

/** Type of cell structure in use */
CellStructure cell_structure;

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

#ifdef ADDITIONAL_CHECKS
static void check_cells_consistency()
{
  int c, index;
  IntList used;
  alloc_intlist(&used, n_cells);
  memset(used.e, 0, n_cells*sizeof(int));
  
  for (c = 0; c < local_cells.n; c++) {
    index = (void *)local_cells.cell[c] - (void *)cells;
    if ((index % sizeof(Cell)) != 0) {
      fprintf(stderr, "%d: local cell pointer not even aligned, certainly wrong (local_cell[%d], index=%d).\n", this_node, c, index);
      errexit();
    }
    index /= sizeof(Cell);
    if (index < 0 || index >= n_cells) {
      fprintf(stderr, "%d: local cell pointer out of range, maybe old leftover (local_cell[%d]).\n", this_node, c);
      errexit();
    }
    if (used.e[index]) {
      fprintf(stderr, "%d: local cell is already pointed to (local_cell[%d]).\n", this_node, c);
      errexit();
    }
    used.e[index] = 1;
  }

  for (c = 0; c < ghost_cells.n; c++) {
    index = (void *)ghost_cells.cell[c] - (void *)cells;
    if ((index % sizeof(Cell)) != 0) {
      fprintf(stderr, "%d: ghost cell pointer not even aligned, certainly wrong (ghost_cell[%d], index=%d).\n", this_node, c, index);
      errexit();
    }
    index /= sizeof(Cell);
    if (index < 0 || index >= n_cells) {
      fprintf(stderr, "%d: ghost cell pointer out of range, maybe old leftover (ghost_cell[%d]).\n", this_node, c);
      errexit();
    }
    if (used.e[index]) {
      fprintf(stderr, "%d: ghost cell is already pointed to (ghost_cell[%d]).\n", this_node, c);
      errexit();
    }
    used.e[index] = 1;
  }
  for (c = 0; c < n_cells; c++)
    if (!used.e[c]) {
      fprintf(stderr, "%d: cell %d is not used anywhere.\n", this_node, c);
      errexit();
    }
  realloc_intlist(&used, 0);
}
#endif

/*@}*/

/************************************************************/
void cells_pre_init()
{
  CellPList tmp_local;
  CELL_TRACE(fprintf(stderr, "%d: cells_pre_init\n",this_node));
  /* her local_cells has to be a NULL pointer */
  if(local_cells.cell != NULL) {
    fprintf(stderr,"Wrong usage of cells_pre_init!\n");
    errexit();
  }
  memcpy(&tmp_local,&local_cells,sizeof(CellPList));
  dd_topology_init(&tmp_local);
}

void realloc_cells(int size)
{
  int i;
  CELL_TRACE(fprintf(stderr, "%d: realloc_cells %d\n", this_node, size));
  /* free all memory associated with cells to be deleted. */
  for(i=size; i<n_cells; i++) {
    realloc_particlelist(&cells[i],0);
  }
  /* resize the cell list */
  if(size != n_cells) {
    cells = (Cell *) realloc(cells, sizeof(Cell)*size);
  }
  /* initialize new cells */
  for(i=n_cells; i<size; i++) {
    init_particleList(&cells[i]);
  }
  n_cells = size;
}  

/************************************************************/
static void topology_release(int cs) {
  switch (cs) {
  case CELL_STRUCTURE_CURRENT:
    topology_release(cell_structure.type);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_release();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_release();
    break;
  default:
    fprintf(stderr, "ERROR: attempting to sort the particles in an unknown way\n");
    errexit();
  }
}

/************************************************************/
static void topology_init(int cs, CellPList *local) {
  switch (cs) {
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type, local);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(local);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init(local);
    break;
  default:
    fprintf(stderr, "ERROR: attempting to sort the particles in an unknown way\n");
    errexit();
  }
}

/************************************************************/
void cells_re_init(int new_cs)
{
  CellPList tmp_local;
  Cell *tmp_cells;
  int tmp_n_cells,i;

  CELL_TRACE(fprintf(stderr, "%d: cells_re_init: convert type (%d->%d)\n",
		     this_node, cell_structure.type, new_cs));

  invalidate_ghosts();

  CELL_TRACE({
    int p;
    for (p = 0; p < n_total_particles; p++)
      if (local_particles[p])
	fprintf(stderr, "%d: cells_re_init: got particle %d\n", this_node, p);
  }
	     );
    
  topology_release(cell_structure.type);
  /* MOVE old local_cell list to temporary buffer */
  memcpy(&tmp_local,&local_cells,sizeof(CellPList));
  init_cellplist(&local_cells);

  /* MOVE old cells to temporary buffer */
  tmp_cells   = cells;
  tmp_n_cells = n_cells;
  cells   = NULL;
  n_cells = 0;

  topology_init(new_cs, &tmp_local);

  /* finally deallocate the old cells */
  realloc_cellplist(&tmp_local,0);
  for(i=0;i<tmp_n_cells;i++)
    realloc_particlelist(&tmp_cells[i],0);

  free(tmp_cells);
  CELL_TRACE(fprintf(stderr, "%d: old cells deallocated\n",this_node));

  CELL_TRACE({
    int p;
    for (p = 0; p < n_total_particles; p++)
      if (local_particles[p])
	fprintf(stderr, "%d: cells_re_init: now got particle %d\n", this_node, p);
  }
	     );

#ifdef ADDITIONAL_CHECKS
  check_cells_consistency();
#endif
}

/*************************************************/
int cells_get_n_particles()
{
  int c, cnt = 0;
  for (c = 0; c < local_cells.n; c++)
    cnt += local_cells.cell[c]->n;
  return cnt;
}

/*************************************************/
int cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv)
{
  if (argc < 1) {
    Tcl_AppendResult(interp, "usage: cellsystem <system> <params>", (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG1_IS_S("domain_decomposition"))
    mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC);
  else if (ARG1_IS_S("nsquare"))
    mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE);
  else {
    Tcl_AppendResult(interp, "unkown cell structure type \"", argv[0],"\"", (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

/*************************************************/
void print_local_particle_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0 ; i < np; i++) {
      fprintf(stderr,"%d: local cell %d contains part id=%d pos=(%f,%f,%f)\n",
	      this_node, c, part[i].p.identity,
	      part[i].r.p[0], part[i].r.p[1], part[i].r.p[2]);
      cnt++;
    }
  }
  fprintf(stderr,"%d: found %d particles\n",this_node,cnt);
}

/*************************************************/
void print_ghost_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0 ; i < np; i++) {
      fprintf(stderr,"%d: local cell %d contains ghost id=%d pos=(%f,%f,%f)\n",
	      this_node, c, part[i].p.identity,
	      part[i].r.p[0], part[i].r.p[1], part[i].r.p[2]);
      cnt++;
    }
  }
  fprintf(stderr,"%d: found %d ghosts\n",this_node,cnt);
}
