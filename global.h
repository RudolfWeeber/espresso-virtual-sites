#ifndef GLOBAL_H
#define GLOBAL_H
#include <tcl.h>

/**********************************************
 * global variables
 * put everything here that should be available
 * to more than one module, especially those
 * variables that should be treated as Tcl
 * datafields.
 * Please mark in which fiel the variable is
 * defined.
 **********************************************/

/****************************************
 * mpi related stuff from communication.c
 *****************************************/
extern int this_node;
extern int nprocs;

/****************************************
 * topology from grid.c
 ****************************************/
extern int processor_grid[3];
extern int pe_pos[3];
extern int neighbors[6];

/****************************************
 * box dimensions from global.c
 *****************************************/
extern double box_l[3];
/* size of local box. */
extern double local_box_l[3];
/* left corner of local box. */
extern double my_left[3];
/* right corner of local box. */
extern double my_right[3];

/****************************************
 * particle data from global.c
 ****************************************/

/** Field to hold particle information
 *  of local particles. */
typedef struct {
  int    identity;
  int    type;

  /* periodically folded position */
  double p[3];
  /* index of the simulation box image where the particle really sits */
  int    i[3];
  double q;

  double v[3];
  double f[3];

  int   n_bonds;
  int max_bonds;
  int    *bonds;
} Particle;

/** size of local particle array. */
extern int   max_particles;
/** number of particles belonging to that node. */
extern int     n_particles;
/** number of ghost particle belonging to that node. */
extern int     n_ghosts;
/** local particle array. */
extern Particle *particles;

/* total number of particles in the system. */
extern int n_total_particles;
/* used only on master node: particle->node mapping */
extern int  *particle_node;

/** Mapping between particle identity and local index. 
 *    You find the local index of particle i 
 *    at position i of this field. 
 *    A particle that is not in the processors domain 
 *    (including its ghostshell) is marked with -1.
 */
extern int *local_index;

/** allocate storage for local particles and ghosts.
    Given size is rounded up to multiples of
    PART_INCREMENT */
void realloc_particles(int size);

/** search for a specific particle, returns field index */
int got_particle(int part);

/** add a particle, returns field index */
int add_particle(int part);

/** allocate space for a particle, returns field index */
int alloc_particle();

/** fold particle coordinates to primary simulation box */
void fold_particle(double pos[3],int image_box[3]);

/** unfold particle coordinates to physical position */
void unfold_particle(double pos[3],int image_box[3]);

/** free a particle */
void free_particle(int index);

/** add particle to particle->node map */
void map_particle_node(int part, int node);

/** rebuild particle->node map from scratch */
void build_particle_node();

/** update n_total_particles on slave nodes
 *  and invalidate particle_node */
void particle_finalize_data();

/****************************************
 * nonbonded interactions from global.c
 ****************************************/

/** number of particle types. */
extern int n_particle_types;

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;

  /* don't know which else, since electrostatic is different...
     but put rest here too. */
} IA_parameters;

/* size is n_particle_types^2 */
extern IA_parameters *ia_params;

/** number of interaction types. */
extern int n_interaction_types;

/** get interaction particles between particle sorts i and j */
IA_parameters *get_ia_param(int i, int j);

/** get interaction particles between particle sorts i and j.
    returns NULL if i or j < 0, allocates if necessary */
IA_parameters *safe_get_ia_param(int i, int j);

/** realloc n_particle_types */
void realloc_ia_params(int nsize);

/** initialize interaction parameters */
void initialize_ia_params(IA_parameters *params);

/** copy interaction parameters */
void copy_ia_params(IA_parameters *dst, IA_parameters *src);

/****************************************
 * bonded interactions from global.c
 ****************************************/

/** possible values for bonded_ia_type */
#define BONDED_IA_DIHEDRAL 0
#define BONDED_IA_ANGLE    1

typedef union {
  int bonded_ia_type;
  struct {
    int dummy;
  } dihedral;
  struct {
    int dummy;
  } angle;
} Bonded_ia_parameters;

/** field defining the bonded ia types */
extern int n_bonded_ia;
extern Bonded_ia_parameters *bonded_ia_params;

/** reallocate particles bonds */
void realloc_bonds(int index, int size);

/****************************************
 * integration from integrator.c
 ****************************************/

/** time step for integration */
extern double time_step;
/** maximal interaction cutoff. */
extern double max_cut;
/** verlet list skin. */
extern double skin;
/** maximal interaction range (max_cut + skin). */
extern double max_range;
/** maximal interaction range squared. */
extern double max_range2;

/** Flag for integrator.
    Wether to calculate the forces before the first step. */
extern int calc_forces_first;

/****************************************
 * Verlet list from verlet.c
 ****************************************/
extern int   n_verletList;
extern int max_verletList;
extern int    *verletList;

/** Flag for rebuilding the verlet list. */
extern int rebuild_verletlist;

/**********************************************
 * description of global variables
 * add any variable that should be handled
 * automatically in global.c. This includes
 * distribution to other nodes and
 * read/user-defined access from Tcl.
 **********************************************/

/** possible field types */
#define TYPE_INT    0
#define TYPE_DOUBLE 1

/** maximal dimension of a writable datafield */
#define MAX_DIMENSION 64

/* set callback procedure */
typedef int (SetCallback)(Tcl_Interp *interp, void *data);

/* variable descriptor */
typedef struct {
  void        *data;      /* physical address */
  int          type;      /* int or double */
  int          dimension; /* field dimension */
  const char  *name;      /* name assigned in Tcl */
  SetCallback *changeproc;/* procedure called if value should be
			     changed. Maybe ro_callback for
			     non-writeable variables */
} Datafield;

extern const Datafield fields[];

/**********************************************
 * misc procedures
 **********************************************/

/** initialize data fields */
void init_data();

/** call if topology (grid, box dim, ...) changed */
void changed_topology();

#endif
