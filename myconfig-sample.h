/** \file myconfig-sample.h

    This is a sample for the file myconfig.h.

    Uncomment any of the following lines to myconfig.h to activate
    the corresponding feature of Espresso. It is recommended to turn
    only those features on that you actually need to optimize the
    performance of Espresso for your problem. For details on these
    features see the user's guide, Sec. 2.6 "Configuration options".

    To access the information on the compilation status of the code
    you are working with in your Espresso Tcl-script, use the
    corresponding \ref tcl_features "Tcl-commands".

    If you add a new feature to Espresso, you also have to add the
    corresponding lines in the function \ref compilation_callback and
    to add documentation in <tt>doc/text/features.doc</tt>.
 */

/**********************************************************************/
/*                       general core features                        */
/**********************************************************************/

/* #define PARTIAL_PERIODIC */
/* #define ELECTROSTATICS */
/* #define ROTATION */
/* #define ROTATIONAL_INTERIA */
/* #define DIPOLES */
/* #define MAGNETOSTATICS */
/* #define EXTERNAL_FORCES */
/* #define CONSTRAINTS */
/* #define MASS */
/* #define EXCLUSIONS */
/* #define COMFORCE */
/* #define COMFIXED */
/* #define MOLFORCES */
/* #define BOND_CONSTRAINT */
/* #define MODES */
/* #define BOND_VIRTUAL */

// To use address, also activate VIRTUAL_SITES and VIRTUAL_SITES_COM
/* #define ADRESS*/


/* #define DAWAANR */
/* #define MAGNETIC_DIPOLAR_DIRECT_SUM */
/* #define MDLC */
/* #define METADYNAMICS */
/* #define OVERLAPPED */

// When using virtual_sites, activate ONE implementation, below
/* #define VIRTUAL_SITES */
// This implementation puts the virtual particle in the center of mass
// of all real particles in the same molecule
/* #define VIRTUAL_SITES_COM */
// This implementation puts a virtual particle relative to a real particle's
// positino and orientation 
/* #define VIRTUAL_SITES_RELATIVE */

// Activate the following, if you don't want the velocity of virtual sites to be calculated
// #define VIRTUAL_SITES_NO_VELOCITY

// Activate the following, if you want the langevin thermostat to act on virtual sites
// #define VIRTUAL_SITES_THERMOSTAT

// Activate the following, if you do NOT want the Langevin thermostat to act on non-virtual particles
// #define THERMOSTAT_IGNORE_NON_VIRTUAL


/**********************************************************************/
/*                        integrator features                         */
/**********************************************************************/

/* #define NEMD */
/* #define NPT */ 
/* #define DPD */
/* #define TRANS_DPD */
/* #define INTER_DPD */

/*Note: Activate ONLY ONE dpd mass  out of the following! */
/* #define DPD_MASS_RED */
/* #define DPD_MASS_LIN */

/* #define LB */

/**********************************************************************/
/*                           interactions                             */
/**********************************************************************/

/* #define TABULATED */
/* #define LENNARD_JONES */
/* #define LJ_WARN_WHEN_CLOSE */
/* #define LENNARD_JONES_GENERIC */
/* #define LJCOS */
/* #define LJCOS2 */
/* #define LJ_ANGLE */
/* #define GAY_BERNE */
/* #define SMOOTH_STEP */
/* #define HERTZIAN */
/* #define BMHTF_NACL */
/* #define MORSE */
/* #define BUCKINGHAM */
/* #define SOFT_SPHERE */
/* #define INTER_DPD */
/* #define INTER_RF */
/* #define MOL_CUT */
/* #define TUNABLE_SLIP */
/* #define NO_INTRA_NB */

/* Note: Activate ONLY ONE bonded angle potential out of the following! */
/* #define BOND_ANGLE_HARMONIC */
/* #define BOND_ANGLE_COSINE */
/* #define BOND_ANGLE_COSSQUARE */

/* #define BOND_ANGLEDIST */
/* #define BOND_ENDANGLEDIST */

/* activate the old dihedral form. Only uncomment this if you need
   to run old code using the previous phase definition. */
/* #define OLD_DIHEDRAL */

/**********************************************************************/
/*                            debugging                               */
/**********************************************************************/

/* #define ADDITIONAL_CHECKS */

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
/* #define PTENSOR_DEBUG*/
/* #define MEM_DEBUG */
/* #define MAGGS_DEBUG */
/* #define LB_DEBUG */
/* #define VIRTUAL_SITES_DEBUG */

/* #define ASYNC_BARRIER */

/* #define MPI_CORE */
/* #define FORCE_CORE */

/* #define OLD_RW_VERSION*/
