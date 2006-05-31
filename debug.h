// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file debug.h
 This file controls debug facilities. 

 <b>Responsible:</b>
 <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

 The implementation is found in
 \ref debug.c "debug.c".
 This file contains a couple of commented out defines. If they are commented in,
 they activate debug output for various aspects:
 <ul>
 <li> \verbatim #define COMM_DEBUG \endverbatim activate printing of communicated actions.
 <li> \verbatim #define INTEG_DEBUG \endverbatim activate integration debug output.
 <li> \verbatim #define CELL_DEBUG \endverbatim activate cell code debug output.
 <li> \verbatim #define GHOST_DEBUG \endverbatim activate ghost code debug output.
 <li> \verbatim #define LATTICE_DEBUG \endverbatim activate lattice code debug output.
 <li> \verbatim #define HALO_DEBUG \endverbatim activate halo code debug output.
 <li> \verbatim #define GRID_DEBUG \endverbatim activate grid debug output.
 <li> \verbatim #define VERLET_DEBUG \endverbatim activate verlet debug output.
 <li> \verbatim #define PARTICLE_DEBUG \endverbatim activate particle data related debug output.
 <li> \verbatim #define P3M_DEBUG \endverbatim activate P3M related debug output.
 <li> \verbatim #define FFT_DEBUG \endverbatim activate FFT related debug output.
 <li> \verbatim #define RANDOM_DEBUG \endverbatim activate debug output for the random number generator.
 <li> \verbatim #define FORCE_DEBUG \endverbatim activate detailed Force debugging -> gives explicitely all forces.
 <li> \verbatim #define THERMO_DEBUG \endverbatim activate debug output for the theromostat.
 <li> \verbatim #define LJ_DEBUG \endverbatim activate debug output for the LJ-potential.
 <li> \verbatim #define MORSE_DEBUG \endverbatim activate debug output for the MORSE-potential.
 <li> \verbatim #define ESR_DEBUG \endverbatim activate ??? debug output.
 <li> \verbatim #define ESK_DEBUG \endverbatim activate ??? debug output.
 <li> \verbatim #define FENE_DEBUG \endverbatim activate debug output for the FENE-potential.
 <li> \verbatim #define GHOST_FORCE_DEBUG \endverbatim activate debug output for ghost forces.
 <li> \verbatim #define ONE_PART_DEBUG <part_id> \endverbatim activate debug output, tracing everything done with particle "<part_id>".
 <li> \verbatim #define STAT_DEBUG \endverbatim activate statistics related debug output.
 <li> \verbatim #define POLY_DEBUG \endverbatim activate polymer setup debug output.
 <li> \verbatim #define MOLFORCES_DEBUG \endverbatim activate molforces debug output.
 <li> \verbatim #define MEM_DEBUG \endverbatim activate memory allocation/freeing logging. WARNING: This generates TONS of
 output. This can only be used reasonably for a few time steps with a few particles.
 <li> \verbatim #define MAGGS_DEBUG \endverbatim activate debugging output of the Maggs electrostatics/dynamics method.
 <li> \verbatim #define LB_DEBUG \endverbatim activate debugging output for Lattice Boltzmann fluid dynamics.
 <li> \verbatim #define MPI_CORE \endverbatim generate a core dump when exiting abnormally due
 to MPI errors.
 <li> \verbatim #define FORCE_CORE \endverbatim generate a core dump even on regular termination.
 <li> \verbatim #define ADDITIONAL_CHECKS \endverbatim adds additional checks which indicate failure
 of the code.
 </ul>

 For every define there exists a macro that can be used to encapsulate short lines (like printf("...",...);)
 of code that should be executed iff the respective *_DEBUG macro is defined.
*/

#include <tcl.h>

/* #define COMM_DEBUG */
/* #define EVENT_DEBUG */
/* #define INTEG_DEBUG */
/* #define CELL_DEBUG */
/* #define GHOST_DEBUG */
/* #define LATTICE_DEBUG */
/* #define HALO_DEBUG */
/* #define GRID_DEBUG */
/* #define VERLET_DEBUG */
/* #define PARTICLE_DEBUG */
/* #define P3M_DEBUG */
/* #define EWALD_DEBUG */
/* #define FFT_DEBUG */
/* #define RANDOM_DEBUG */
/* #define FORCE_DEBUG */
/* #define THERMO_DEBUG */ 
/* #define LJ_DEBUG */
/* #define MORSE_DEBUG */
/* #define ESR_DEBUG */
/* #define ESK_DEBUG */
/* #define FENE_DEBUG */
/* #define GHOST_FORCE_DEBUG */
/* #define ONEPART_DEBUG 13 */
/* #define STAT_DEBUG */ 
/* #define POLY_DEBUG */
/* #define MOLFORCES_DEBUG */
/* #define MEM_DEBUG */
/* #define MAGGS_DEBUG */
/* #define LB_DEBUG */

/** do an MPI_Barrier at the beginning of the communication slave loop.
    Helps a lot in debugging */
/* #define ASYNC_BARRIER */

/** generate a core file on MPI errors */
#define MPI_CORE
/** generate a core file on errexit() */
#define FORCE_CORE

/* #define ADDITIONAL_CHECKS */

#ifdef MEM_DEBUG
#ifdef __GNUC__
#define realloc(v,s) __realloc((v),(s),__FILE__, __LINE__)
#define malloc(s) __malloc((s),__FILE__, __LINE__)
#define free(v) __free((v),__FILE__, __LINE__)
#else
#define realloc(v,s) __realloc((v),(s), "no line info", 0)
#define malloc(s) __malloc((s), "no line info", 0)
#define free(v) __free((v),"no line info", 0)
#endif

/** memory allocation test routine */
void *__realloc(void *old, unsigned int size, char *where, int line);

/** memory allocation test routine */
void *__malloc(unsigned int size, char *where, int line);

/** memory allocation test routine */
void __free(void *p, char *where, int line);

#endif

/** callback for debug status. */
int debug_callback(Tcl_Interp *interp);

#if defined FORCE_CORE || defined MPI_CORE
/** this functions kills the task with SIGSEGV */
void core();
#endif

#ifdef ADDITIONAL_CHECKS

/** this performs a lot of tests which will very likely detect corruptions of
    \ref local_particles and the cell structure.
*/
void check_particle_consistency();

/** check the consistency of the cells and particle_node. Called from
    mpi_bcast_event(CHECK_PARTICLES)
*/
void check_particles();

#endif

/** Print all particle positions contained in \ref cells::cells array. */
void print_particle_positions();
/** Print all particle forces contained in \ref cells::cells array. */
void print_particle_forces();

/** by setting this variable to 1, a regular exit is
    indicated. In that case, no core dump is generated.
*/
extern int regular_exit;

/** Identity of the particle to chekc extenively if ONEPART_DEBUG is defined. */
extern int check_id;

#ifdef COMM_DEBUG
#define COMM_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff COMM_DEBUG is set. */
#define COMM_TRACE(cmd)
#endif

#ifdef EVENT_DEBUG
#define EVENT_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff EVENT_DEBUG is set. */
#define EVENT_TRACE(cmd)
#endif

#ifdef PARTICLE_DEBUG
#define PART_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff PARTICLE_DEBUG is set. */
#define PART_TRACE(cmd)
#endif

#ifdef INTEG_DEBUG
#define INTEG_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff INTEG_DEBUG is set. */
#define INTEG_TRACE(cmd)
#endif

#ifdef CELL_DEBUG
#define CELL_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff CELL_DEBUG is set. */
#define CELL_TRACE(cmd)
#endif

#ifdef GHOST_DEBUG
#define GHOST_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff GHOST_DEBUG is set. */
#define GHOST_TRACE(cmd)
#endif

#ifdef HALO_DEBUG
#define HALO_TRACE(cmd) { cmd; }
#else
/** Equals { cmd  } iff HALO_DEBUG is set. */
#define HALO_TRACE(cmd)
#endif

#ifdef GRID_DEBUG
#define GRID_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff GRID_DEBUG is set. */
#define GRID_TRACE(cmd)
#endif

#ifdef LATTICE_DEBUG
#define LATTICE_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff LATTICE_DEBUG is set. */
#define LATTICE_TRACE(cmd)
#endif

#ifdef FORCE_DEBUG
#define FORCE_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff FORCE_DEBUG is set. */
#define FORCE_TRACE(cmd)
#endif

#ifdef VERLET_DEBUG
#define VERLET_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff VERLET_DEBUG is set. */
#define VERLET_TRACE(cmd)
#endif

#ifdef P3M_DEBUG
#define P3M_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff P3M_DEBUG is set. */
#define P3M_TRACE(cmd)
#endif

#ifdef EWALD_DEBUG
#define EWALD_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } if EWALD_DEBUG is set. */
#define EWALD_TRACE(cmd)
#endif

#ifdef MAGGS_DEBUG
#define MAGGS_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff MAGGS_DEBUG is set. */
#define MAGGS_TRACE(cmd)
#endif

#ifdef FFT_DEBUG
#define FFT_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff FFT_DEBUG is set. */
#define FFT_TRACE(cmd)
#endif

#ifdef RANDOM_DEBUG
#define RANDOM_TRACE(cmd) { cmd; }
#else
#define RANDOM_TRACE(cmd)
#endif

#ifdef THERMO_DEBUG
#define THERMO_TRACE(cmd) { cmd; }
#else
#define THERMO_TRACE(cmd)
#endif

#ifdef LJ_DEBUG
#define LJ_TRACE(cmd) { cmd; }
#else
#define LJ_TRACE(cmd)
#endif

#ifdef MORSE_DEBUG
#define MORSE_TRACE(cmd) { cmd; }
#else
#define MORSE_TRACE(cmd)
#endif

#ifdef BUCK_DEBUG
#define BUCK_TRACE(cmd) { cmd; }
#else
#define BUCK_TRACE(cmd)
#endif

#ifdef ESR_DEBUG
#define ESR_TRACE(cmd) { cmd; }
#else
#define ESR_TRACE(cmd)
#endif

#ifdef ESK_DEBUG
#define ESK_TRACE(cmd) { cmd; }
#else
#define ESK_TRACE(cmd)
#endif

#ifdef FENE_DEBUG
#define FENE_TRACE(cmd) { cmd; }
#else
#define FENE_TRACE(cmd)
#endif

#ifdef GHOST_FORCE_DEBUG
#define GHOST_FORCE_TRACE(cmd) { cmd; }
#else
#define GHOST_FORCE_TRACE(cmd)
#endif

#ifdef ONEPART_DEBUG
#define ONEPART_TRACE(cmd) { cmd; }
#else
#define ONEPART_TRACE(cmd)
#endif

#ifdef STAT_DEBUG
#define STAT_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff STAT_DEBUG is set. */
#define STAT_TRACE(cmd)
#endif

#ifdef POLY_DEBUG
#define POLY_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff POLY_DEBUG is set. */
#define POLY_TRACE(cmd)
#endif

#ifdef MOLFORCES_DEBUG
#define MOLFORCES_TRACE(cmd) { cmd; }
#else
#define MOLFORCES_TRACE(cmd)
#endif

#ifdef LB_DEBUG
#define LB_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff LB_DEBUG is set. */
#define LB_TRACE(cmd)
#endif
