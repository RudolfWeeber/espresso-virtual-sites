// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>

*/

#include <tcl.h>
#include "particle_data.h"
#include "interaction_data.h"
#include "utils.h"
#include "topology.h"

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {
  /** Status flag for observable calculation. 
      For 'analyze energy': 0 re-initialize observable struct, else every thing is fine, calculation can start. 
      For 'analyze pressure' and 'analyze p_inst': 0 or !(1+v_comp) re-initialize, else all OK. */
  int init_status;

  /** Array for observables on each node. */
  DoubleList data;

  /** number of coulomb interactions */
  int n_coulomb;
  /** number of non bonded interactions */
  int n_non_bonded;

  /** start of bonded interactions. Right after the special ones */
  double *bonded;
  /** start of observables for non-bonded interactions. */
  double *non_bonded;
  /** start of observables for coulomb interaction. */
  double *coulomb;

  /** number of doubles per data item */
  int chunk_size;
} Observable_stat;

/** Structure used only in the pressure and stress tensor calculation to distinguish 
    non-bonded intra- and inter- molecular contributions. */
typedef struct {
  /** Status flag for observable calculation.
      For 'analyze energy': 0 re-initialize observable struct, else every thing is fine, calculation can start.
      For 'analyze pressure' and 'analyze p_inst': 0 or !(1+v_comp) re-initialize, else all OK. */
  int init_status_nb;

  /** Array for observables on each node. */
  DoubleList data_nb;

  /** number of non bonded interactions */
  int n_nonbonded;

  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

  /** number of doubles per data item */
  int chunk_size_nb;
} Observable_stat_non_bonded;

/*@}*/

/** \name Exported Variables
    Previous particle configurations (needed for offline analysis
    and correlation analysis in \ref #analyze)
*/
/************************************************************/
/*@{*/
extern double **configs;
extern int n_configs;
extern int n_part_conf;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Implements the Tcl command \ref tcl_analyze. This allows for basic system analysis,
    both online and offline.
*/
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** the minimal distance of two particles with types in set1 rsp. set2.
    @param set1 types of particles
    @param set2 types of particles
    @return the minimal distance of two particles */
double mindist(IntList *set1, IntList *set2);

/** calculate the aggregate distribution for molecules.
    @param dist_criteria2 distance criteria squared
    @param min_contact minimum number of contacts 
    @param s_mol_id start molecule id
    @param f_mol_id finish molecule id
    @param head_list 
    @param link_list
    @param agg_id_list
    @param agg_num
    @param agg_size
    @param agg_max
    @param agg_min
    @param agg_avg
    @param agg_std
    @param charge_criteria
*/
int aggregation(double dist_criteria2, int min_contact, int s_mol_id, int f_mol_id, int *head_list, int *link_list,
		int *agg_id_list, int *agg_num, int *agg_size, int *agg_max, int *agg_min, int *agg_avg, int *agg_std, int charge_criteria);

/** returns all particles within a given radius r_catch around a position.
    @param pos position of sphere of point
    @param r_catch the radius around the position
    @param il the list where to store the particles indices */
void nbhood(double pos[3], double r_catch, IntList *il, int planedims[3]);

/** minimal distance to point.
    @param pos point
    @param pid  if a valid particle id, this particle is omitted from minimization
                (this is a good idea if the posx, posy, posz is the position of a particle).
    @return the minimal distance of a particle to coordinates (\<posx\>, \<posy\>, \<posz\>). */
double distto(double pos[3], int pid);

/** numerical solution for the integration constant \f$\gamma\f$ in the cell model, determined by 
    \f[\gamma\,\ln\frac{R}{r_0}=\arctan\frac{1}{\gamma}+\arctan\frac{\xi_M-1}{\gamma}\f]
    from which the second integration constant, the Manning radius \f$R_M\f$, follows to
    \f[R_M = R\cdot\exp\left(-\frac{1}{\gamma}\cdot\arctan\frac{1}{\gamma}\right)\f]
    Any value \f$\xi_M\geq 0\f$ is allowed, the function will automatically ensure the 
    analytical continuation required for \f$\xi_M<\ln(R/r_0)/(1+\ln(R/r_0))\f$, in which case 
    \f$\gamma\f$ becomes imaginary.
    @param xi_m   Manning parameter \f$\xi_M=\ell_B/a\f$ (with Bjerrum-length \f$\ell_B\f$ and charge distance \f$a\f$)
    @param Rc     outer radius \f$R_C\f$ of the cylindrical cell around each polyelectrolyte
    @param ro     inner radius \f$r_0\f$ of the cylindrical cell around each polyelectrolyte
    @param gacc   the accuracy up to which \f$\gamma\f$ should be determined
    @param maxtry maximum number of interations to find a solution 
    @param result pointer to double array containing \f$\gamma\f$ and \f$R_M\f$, 
                  and a third entry which is -1.0 if \f$\gamma\f$ is imaginary, +1.0 else. */
void calc_cell_gpb(double xi_m, double Rc, double ro, double gacc, int maxtry, double *result);

/** appends particles' positions in 'partCfg' to \ref #configs */
void analyze_append();

/** appends the configuration stored in 'config[3*count]' to \ref #configs
    @param config the configuration which should be added 
    @param count  how many particles in 'config' */
void analyze_configs(double *config, int count);

/** removes configs[0], pushes all entries forward, appends current 'partCfg' to last spot */
void analyze_push();

/** replaces configs[ind] with current 'partCfg'
    @param ind the entry in \ref #configs to be replaced */
void analyze_replace(int ind);

/** removes configs[ind] and shrinks the array accordingly
    @param ind the entry in \ref #configs to be removed */
void analyze_remove(int ind);

/** Calculates the distribution of particles around others. 
    Calculates the distance distribution of particles with types given
    in the p1_types list around particles with types given in the
    p2_types list. The distances range from r_min to r_max, binned
    into r_bins bins which are either aequidistant (log_flag==0)or
    logarithmically aequidistant (log_flag==1). The result is stored
    in the array dist.
    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param log_flag Wether the bins are (logarithmically) aequidistant.
    @param low      particles closer than r_min
    @param dist     Array to store the result (size: r_bins).
 */
void calc_part_distribution(int *p1_types, int n_p1, int *p2_types, int n_p2, 
			    double r_min, double r_max, int r_bins, int log_flag,
			    double *low, double *dist);
/** Calculates the radial distribution function.

    Calculates the radial distribution function of particles with
    types given in the p1_types list around particles with types given
    in the p2_types list. The range is given by r_min and r_max and
    the distribution function is binned into r_bin bins, which are
    equidistant. The result is stored in the array rdf.

    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param rdf     Array to store the result (size: r_bins).
*/
void calc_rdf(int *p1_types, int n_p1, int *p2_types, int n_p2, 
	      double r_min, double r_max, int r_bins, double *rdf);


/** Calculates the radial distribution function averaged over last n_conf configurations.

    Calculates the radial distribution function of particles with
    types given in the p1_types list around particles with types given
    in the p2_types list. The range is given by r_min and r_max and
    the distribution function is binned into r_bin bins, which are
    equidistant. The result is stored in the array rdf.

    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param rdf      Array to store the result (size: r_bins).
    @param n_conf   Number of configurations from the last stored configuration.
*/
void calc_rdf_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
	      double r_min, double r_max, int r_bins, double *rdf, int n_conf);

/** Calculates the intermolecular radial distribution function averaged over last n_conf configurations.

    Calculates the radial distribution function of particles with
    types given in the p1_types list around particles with types given
    in the p2_types list. The range is given by r_min and r_max and
    the distribution function is binned into r_bin bins, which are
    equidistant. The result is stored in the array rdf.

    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param rdf      Array to store the result (size: r_bins).
    @param n_conf   Number of configurations from the last stored configuration.
*/

void calc_rdf_intermol_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
	      double r_min, double r_max, int r_bins, double *rdf, int n_conf);

/** Calculates the spherically averaged structure factor.

    Calculates the spherically averaged structure factor of particles of a
    given type. The possible wave vectors are given by q = 2PI/L sqrt(nx^2 + ny^2 + nz^2).
    The S(q) is calculated up to a given length measured in 2PI/L (the recommended order of
    the wave vector is less than 20)
    
    @param type   the type of the particles to be analyzed
    @param order  the maximum wave vector length in 2PI/L
    @param sf     array to store the result (size: order^2+1).
*/

void calc_structurefactor(int type, int order, double *sf);
	      
/** returns the minimal squared distance between two positions in the perhaps periodic
    simulation box.
 *  \param pos1  Position one.
 *  \param pos2  Position two.
 */
double min_distance2(double pos1[3], double pos2[3]);

/** returns the minimal distance between two positions in the perhaps periodic
    simulation box.
 *  \param pos1  Position one.
 *  \param pos2  Position two.
 */
MDINLINE double min_distance(double pos1[3], double pos2[3]) {
  return sqrt(min_distance2(pos1, pos2));
}

void centermass(int type, double *com);
void momentofinertiamatrix(int type, double *MofImatrix);
void calculate_verlet_neighbors();

MDINLINE double *obsstat_bonded(Observable_stat *stat, int j)
{
  return stat->bonded + stat->chunk_size*j;
}

MDINLINE double *obsstat_nonbonded(Observable_stat *stat, int p1, int p2)
{
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded + stat->chunk_size*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

MDINLINE double *obsstat_nonbonded_intra(Observable_stat_non_bonded *stat, int p1, int p2)
{
/*  return stat->non_bonded_intra + stat->chunk_size*1; */
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded_intra + stat->chunk_size_nb*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

MDINLINE double *obsstat_nonbonded_inter(Observable_stat_non_bonded *stat, int p1, int p2)
{
/*  return stat->non_bonded_inter + stat->chunk_size*1; */
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded_inter + stat->chunk_size_nb*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

void invalidate_obs();

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded, int n_non_bonded,
			       int n_coulomb, int chunk_size);

void obsstat_realloc_and_clear_non_bonded(Observable_stat_non_bonded *stat_nb, int n_nonbonded, int chunk_size_nb);

/*@}*/

#endif
