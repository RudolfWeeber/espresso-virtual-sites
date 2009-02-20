// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file p3m-dipoles.c  P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignement.c  "p3m-assignement.c"
 
 NB: In general the magnetic dipole-dipole functions bear the same name than the charge-charge but,
     adding in front of the name a D   and replacing where "charge" appears by "dipole". In this way
     one can recognize the similarity of the functions but avoiding nasty confusions in their use.

 PS: By default the magnetic epsilon is 1, it has no sense to change the epsilon in magnetic cases, 
     but we left along the code the epsilon handling because in this way in the future it will be
     easier to map to the electrical dipoles. Please do not get rid of the epsilon along the code.  
*/

#ifdef MAGNETOSTATICS

/* only include from within p3m.c */
#ifndef P3M_C_CURRENT
#error never compile this file file directly, it is part of p3m.c
#endif
  
/************************************************
 * DEFINES
 ************************************************/

/* MPI tags for the charge-charge p3m communications: */
/** Tag for communication in P3M_init() -> send_calc_mesh(). */
#define REQ_P3M_INIT_D   2001
/** Tag for communication in gather_fft_grid(). */
#define REQ_P3M_GATHER_D 2011
/** Tag for communication in spread_force_grid(). */
#define REQ_P3M_SPREAD_D 2021


/********** definition of Variables *************/
/** interpolation of the charge assignment function. */
  double *Dint_caf[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};

/** position shift for calc. of first assignment mesh point. */
  double Dpos_shift;

/** help variable for calculation of aliasing sums */  
  double *Dmeshift = NULL;
  
 /** Spatial differential operator in k-space. We use an i*k differentiation. */ 
  double *Dd_op = NULL;

  /** Force optimised influence function (k-space) */
  double *Dg_force = NULL;
  
   /** Energy optimised influence function (k-space) */
  double *Dg_energy = NULL;
  
   /** number of magnetic particles on the node. */
  int Dca_num=0;

  /** Charge fractions for mesh assignment. */
  double *Dca_frac = NULL;
  
  /** index of first mesh point for charge assignment. */  
  int *Dca_fmp = NULL;
  
  /** number of permutations in k_space */  
  int Dks_pnum;

  /** Sum of square of magnetic dipoles (only on master node). */
  double p3m_sum_mu2 = 0.0;
  
  /** number of dipolar particles (only on master node). */
  int p3m_sum_dip_part=0; 

  /** local mesh. */
  local_mesh Dlm;
  /** send/recv mesh sizes */
  send_mesh  Dsm;

  /** size of linear array for local CA/FFT mesh . */
  int    Dca_mesh_size;
  
  /** real space mesh (local) for CA/FFT of the dipolar field.*/
  double *Drs_mesh_dip[3] = {NULL,NULL,NULL};

   /** real space mesh (local) for CA/FFT.*/
   double *Drs_mesh = NULL;
  
   /** k space mesh (local) for k space calculation and FFT.*/
   double *Dks_mesh = NULL;

  /** Field to store grid points to send. */
  double *Dsend_grid = NULL; 
  
  /** Field to store grid points to recv */ 
  double *Drecv_grid = NULL;
  
  /** Allocation size of Dsend_grid and Drecv_grid. */
   int Dsend_recv_grid_size=0;
   
   /**Flag to know if we should calculate the constants for the energy 
    (If you neither compute the energy, is a waste of time spendig circa 3 or 4 min computing such constants **/
   int  Dflag_constants_energy_dipolar=0;


/** \name Private Functions */
/************************************************************/
/*@{*/



/** Calculates for magnetic dipoles the properties of the send/recv sub-meshes of the local FFT mesh. 
 *  In order to calculate the recv sub-meshes there is a communication of 
 *  the margins between neighbouring nodes. */ 
void Dcalc_send_mesh();


/** Initializes for magnetic dipoles the (inverse) mesh constant \ref p3m_struct::a (\ref p3m_struct::ai) 
    and the cutoff for charge assignment \ref p3m_struct::cao_cut, which has to be
    done by \ref P3M_init once and by \ref P3M_scaleby_box_l_dipoles whenever the \ref box_l changed.
*/
void DP3M_init_a_ai_cao_cut();


/** Calculate for magnetic dipoles the spacial position of the left down mesh point of the local mesh, to be
    stored in \ref local_mesh::ld_pos; function called by \ref Dcalc_local_ca_mesh once
    and by \ref P3M_scaleby_box_l_dipoles whenever the \ref box_l changed. */
void Dcalc_lm_ld_pos();


/** Gather FFT grid.
 *  After the charge assignment Each node needs to gather the
 *  information for the FFT grid in his spatial domain.
 */
void Dgather_fft_grid(double* mesh);

/** Spread force grid.
 *  After the k-space calculations each node needs to get all force
 *  information to reassigne the forces from the grid to the
 *  particles.
 */
void Dspread_force_grid(double* mesh);

/** realloc charge assignment fields. */
void Drealloc_ca_fields(int newsize);


/** Initializes the (inverse) mesh constant \ref p3m_struct::a (\ref p3m_struct::ai) 
    and the cutoff for charge assignment \ref p3m_struct::cao_cut, which has to be
    done by \ref P3M_init once and by \ref P3M_scaleby_box_l_dipoles whenever the \ref box_l changed.
*/
void DP3M_init_a_ai_cao_cut();


/** checks for correctness for dipoles in P3M of the cao_cut, necessary when the box length changes */
int DP3M_sanity_checks_boxl();


/** Calculate the spacial position of the left down mesh point of the local mesh, to be
    stored in \ref local_mesh::ld_pos; function called by \ref Dcalc_local_ca_mesh once
    and by \ref P3M_scaleby_box_l_dipoles whenever the \ref box_l changed. */
void Dcalc_lm_ld_pos();


/** Calculates properties of the local FFT mesh for the 
    charge assignment process. */
void Dcalc_local_ca_mesh();

/** Interpolates the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
void interpolate_dipole_assignment_function();

/** shifts the mesh points by mesh/2 */
void Dcalc_meshift();

/** Calculates the Fourier transformed differential operator.  
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
void Dcalc_differential_operator();

/** Calculates the influence function optimized for the dipolar forces. */
void Dcalc_influence_function_force();

/** Calculates the influence function optimized for the dipolar energy and torques. */
void Dcalc_influence_function_energy();

/** Calculates the constants necessary to correct the dipolar energy to minimize the error. */
 void compute_constants_energy_dipolar();
 
/** Calculates the aliasing sums for the optimal influence function.
 *
 * Calculates the aliasing sums in the nominator and denominator of
 * the expression for the optimal influence function (see
 * Hockney/Eastwood: 8-22, p. 275).  
 *
 * \param  n           n-vector for which the aliasing sum is to be performed.
 * \param  nominator   aliasing sums in the nominator.
 * \return denominator aliasing sum in the denominator
 */
MDINLINE double Dperform_aliasing_sums_force(int n[3], double nominator[1]);
MDINLINE double Dperform_aliasing_sums_energy(int n[3], double nominator[1]);
/*@}*/




/* Compute the dipolar surface terms */
double calc_surface_term(int force_flag, int energy_flag);


/** \name P3M Tuning Functions (private)*/
/************************************************************/
/*@{*/


// These 3 functions are to tune the P3M code in the case of dipolar interactions

double P3M_DIPOLAR_real_space_error(double box_size, double prefac, double r_cut_iL,
			    int n_c_part, double sum_q2, double alpha_L);
double P3M_DIPOLAR_k_space_error(double box_size, double prefac, int mesh,
			 int cao, int n_c_part, double sum_q2, double alpha_L); 
			 
void P3M_DIPOLAR_tune_aliasing_sums(int nx, int ny, int nz, 
			    int mesh, double mesh_i, int cao, double alpha_L_i, 
			    double *alias1, double *alias2)	;		 

// To compute the value of alpha  through a bibisection method from the formula 33 of JCP115,6351,(2001).
double JJ_rtbisection(double box_size, double prefac, double r_cut_iL, int n_c_part, double sum_q2,  double x1, double x2, double xacc, double tuned_accuracy);


/************************************************************/
/* functions related to the correction of the dipolar p3m-energy */

double  Dipolar_energy_correction;  /* Stores the value of the energy correction due to MS effects */

double JJ_sumi1(double alpha_L);
double JJ_sumi2(double alpha_L);


double P3M_Average_dipolar_SelfEnergy(double box_l, int mesh);
double perform_aliasing_sums_dipolar_self_energy(int n[3]);



/************************************************************/
/* functions related to the correction of the dipolar p3m-energy */


// Do the sum over k<>0 where k={kx,ky,kz} with kx,ky,kz INTEGERS, of
// exp(-PI**2*k**2/alpha**2/L**2)
double JJ_sumi1(double alpha_L){
       int k2,kx,ky,kz,kx2,ky2,limit=60;
       double suma,alpha_L2;
       
       alpha_L2= alpha_L* alpha_L;
       
       //fprintf(stderr,"alpha_L=%le\n",alpha_L); 
       //fprintf(stderr,"PI=%le\n",PI); 
       
       
       suma=0.0;
       for(kx=-limit;kx<=limit;kx++){
         kx2=kx*kx;
       for(ky=-limit;ky<=limit;ky++){
         ky2=ky*ky;
       for(kz=-limit;kz<=limit;kz++){
           k2=kx2+ky2+kz*kz;
           suma+=exp(-PI*PI*k2/(alpha_L*alpha_L));
       }}} 
       suma-=1; //It's easier to substract the term k=0 later than put an if inside the loops
       
       
         //fprintf(stderr,"suma=%le\n",suma); 
     
       
   return suma;
}

/************************************************************/


// Do the sum over n<>0 where n={nx*L,ny*L,nz*L} with nx,ny,nz INTEGERS, of
// exp(-alpha_iL**2*n**2)
double JJ_sumi2(double alpha_L){
       int n2,nx,ny,nz,nx2,ny2,limit=60;
       double suma;
       
 
       
       suma=0.0;
       for(nx=-limit;nx<=limit;nx++){
         nx2=nx*nx;
       for(ny=-limit;ny<=limit;ny++){
         ny2=ny*ny;
       for(nz=-limit;nz<=limit;nz++){
           n2=nx2+ny2+nz*nz;
           suma+=exp(-alpha_L*alpha_L*n2);
       }}} 
       suma-=1; //It's easier to substract the term n=0 later than put an if inside the loops
       
       
       
   return suma;
}




/* ====== Subroutines to compute analyticaly <Uk_p3m> and parse the output .============*/
double P3M_Average_dipolar_SelfEnergy(double box_l, int mesh) {
	int	i,ind,n[3];
	double node_phi = 0.0, phi = 0.0;
	double U2;
	
        int end[3];
        int size=1;
	
	
   for(i=0;i<3;i++) {
    size *= Dfft_plan[3].new_mesh[i];
    end[i] = Dfft_plan[3].start[i] + Dfft_plan[3].new_mesh[i];
  }
 
  
  for(n[0]=Dfft_plan[3].start[0]; n[0]<end[0]; n[0]++){
    for(n[1]=Dfft_plan[3].start[1]; n[1]<end[1]; n[1]++){
      for(n[2]=Dfft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-Dfft_plan[3].start[2]) + Dfft_plan[3].new_mesh[2] *
	((n[1]-Dfft_plan[3].start[1]) + (Dfft_plan[3].new_mesh[1]*(n[0]-Dfft_plan[3].start[0])));

	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	 node_phi += 0.0;
	else if( (n[0]%(p3m.Dmesh[0]/2)==0) &&
		 (n[1]%(p3m.Dmesh[0]/2)==0) &&
		 (n[2]%(p3m.Dmesh[0]/2)==0) )
	  node_phi += 0.0;
	else {
		  U2 = perform_aliasing_sums_dipolar_self_energy(n);
		  node_phi += Dg_energy[ind] * U2*(SQR(Dd_op[n[0]])+SQR(Dd_op[n[1]])+SQR(Dd_op[n[2]]));
	}
      }}}
  
      
     MPI_Reduce(&node_phi, &phi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
     
     phi*=PI/3./box_l/pow(mesh,3);
     

   return phi ;
}



double perform_aliasing_sums_dipolar_self_energy(int n[3])
{
  double u_sum = 0.0;
  /* lots of temporary variables... */
  double f1,sx,sy,sz,mx,my,mz,nmx,nmy,nmz;
  int    limit=P3M_BRILLOUIN+5;

  f1 = 1.0/(double)p3m.Dmesh[0];

  for(mx = -limit; mx <=limit; mx++) {
    nmx = Dmeshift[n[0]] + p3m.Dmesh[0]*mx;
    sx  = pow(sinc(f1*nmx),2.0*p3m.Dcao);
    for(my = -limit; my <= limit; my++) {
      nmy = Dmeshift[n[1]] + p3m.Dmesh[0]*my;
      sy  = sx*pow(sinc(f1*nmy),2.0*p3m.Dcao);
      for(mz = -limit; mz <=limit; mz++) {
	nmz = Dmeshift[n[2]] + p3m.Dmesh[0]*mz;
	sz  = sy*pow(sinc(f1*nmz),2.0*p3m.Dcao);
	u_sum += sz;
      }
    }
  }
  return u_sum;
}




/******************  functions related to the parsing&tuning  of the dipolar parameters **********/
			 

void Dp3m_set_tune_params(double r_cut, int mesh, int cao,
			 double alpha, double accuracy, int n_interpol)
{
  if (r_cut >= 0) {
    p3m.Dr_cut    = r_cut;
    p3m.Dr_cut_iL = r_cut*box_l_i[0];
  }

  if (mesh >= 0)
    p3m.Dmesh[2] = p3m.Dmesh[1] = p3m.Dmesh[0] = mesh;

  if (cao >= 0)
    p3m.Dcao = cao;

  if (alpha >= 0) {
    p3m.Dalpha   = alpha;
    p3m.Dalpha_L = alpha*box_l[0];
  }

  if (accuracy >= 0)
    p3m.Daccuracy = accuracy;

  if (n_interpol != -1)
    p3m.Dinter = n_interpol;

  coulomb.Dprefactor = (temperature > 0) ? temperature*coulomb.Dbjerrum : coulomb.Dbjerrum;

}


/*****************************************************************************/

int Dp3m_set_params(double r_cut, int mesh, int cao,
		   double alpha, double accuracy)
{
  if(r_cut < 0)
    return -1;

  if(mesh < 0)
    return -2;

  if(cao < 1 || cao > 7 || cao > mesh)
    return -3;

  p3m.Dr_cut    = r_cut;
  p3m.Dr_cut_iL = r_cut*box_l_i[0];
  p3m.Dmesh[2]  = p3m.Dmesh[1] = p3m.Dmesh[0] = mesh;
  p3m.Dcao      = cao;

  if (alpha > 0) {
    p3m.Dalpha   = alpha;
    p3m.Dalpha_L = alpha*box_l[0];
  }
  else
    if (alpha != -1.0)
      return -4;

  if (accuracy >= 0)
    p3m.Daccuracy = accuracy;
  else
    if (accuracy != -1.0)
      return -5;

  mpi_bcast_coulomb_params();

  return 0;
}

/*****************************************************************************/


int Dp3m_set_mesh_offset(double x, double y, double z)
{
  if(x < 0.0 || x > 1.0 ||
     y < 0.0 || y > 1.0 ||
     z < 0.0 || z > 1.0 )
    return TCL_ERROR;

  p3m.Dmesh_off[0] = x;
  p3m.Dmesh_off[1] = y;
  p3m.Dmesh_off[2] = z;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}


/*****************************************************************************/
/* We left the handling of the epsilon, duer to portability reasons in the future for the electrical dipoles,
or if people wants to do electrical dipoles alone using the magnetic code .. */

int Dp3m_set_eps(double eps)
{
  p3m.Depsilon = eps;

  fprintf(stderr,">> p3m.Depsilon =%lf\n",p3m.Depsilon);
  fprintf(stderr,"if you are doing true MAGNETIC CALCULATIONS the value of Depsilon should be 1, if you change it, you go on your own risk ...\n");

  mpi_bcast_coulomb_params();

  return TCL_OK;
}


/*****************************************************************************/



int Dp3m_set_ninterpol(int n)
{
  if (n < 0)
    return TCL_ERROR;

  p3m.Dinter = n;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

/*****************************************************************************/




int Dinter_parse_p3m_tune_params(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
  int mesh = -1, cao = -1, n_interpol = -1;
  double r_cut = -1, accuracy = -1;

  while(argc > 0) {
    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut >= -1)) {
	Tcl_AppendResult(interp, "r_cut expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("mesh")) {
      if(! (argc > 1 && ARG1_IS_I(mesh) && mesh >= -1)) {
	Tcl_AppendResult(interp, "mesh expects an integer >= -1",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(cao) && cao >= -1 && cao < 7)) {
	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
			 (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("accuracy")) {
      if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0)) {
	Tcl_AppendResult(interp, "accuracy expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else if (ARG0_IS_S("n_interpol")) {
      if (! (argc > 1 && ARG1_IS_I(n_interpol) && n_interpol >= 0)) {
	Tcl_AppendResult(interp, "n_interpol expects an nonnegative integer",
			 (char *) NULL);
	return TCL_ERROR;
      }
    }
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }
  
  Dp3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy, n_interpol);

  /* check for optional parameters */
  if (argc > 0) {
    if (Dinter_parse_p3m_opt_params(interp, argc, argv) == TCL_ERROR)
      return TCL_ERROR;
  }

  if (adaptive) {
    if(DP3M_adaptive_tune_parameters(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }
  else {
    if(DP3M_tune_parameters(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }

  return TCL_OK;
}



/*****************************************************************************/



int Dinter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh, cao, i;

  if (coulomb.Dmethod != DIPOLAR_P3M && coulomb.Dmethod != DIPOLAR_MDLC_P3M)
    coulomb.Dmethod = DIPOLAR_P3M;
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with dipolar P3M",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (argc < 1) {
    Tcl_AppendResult(interp, "expected: inter dipolar <bjerrum> p3m tune | <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    Tcl_AppendResult(interp, "Node grid not suited for dipolar P3M. Node grid must be sorted, largest first.", (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return Dinter_parse_p3m_tune_params(interp, argc-1, argv+1, 0);

  if (ARG0_IS_S("tunev2"))
    return Dinter_parse_p3m_tune_params(interp, argc-1, argv+1, 1);
      
  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc < 3 || argc > 5) {
    Tcl_AppendResult(interp, "wrong # arguments: inter dipolar <bjerrum> p3m <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if((! ARG_IS_I(1, mesh)) || (! ARG_IS_I(2, cao))) {
    Tcl_AppendResult(interp, "integer expected", (char *) NULL);
    return TCL_ERROR;
  }
	
  if(argc > 3) {
    if(! ARG_IS_D(3, alpha))
      return TCL_ERROR;
  }
  else {
    Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc > 4) {
    if(! ARG_IS_D(4, accuracy)) {
      Tcl_AppendResult(interp, "double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if ((i = Dp3m_set_params(r_cut, mesh, cao, alpha, accuracy)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
      break;
    case -3:
      Tcl_AppendResult(interp, "cao must be between 1 and 7 and less than mesh",
		       (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "accuracy must be positive", (char *) NULL);
      break;
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  }

  return TCL_OK;
}


/*****************************************************************************/



int Dinter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int i; double d1, d2, d3;

  Tcl_ResetResult(interp);

  while (argc > 0) {
    /* p3m parameter: inter */
    if (ARG0_IS_S("n_interpol")) {
      
      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG1_IS_I(i)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 INTEGER parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (Dp3m_set_ninterpol(i) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " argument must be positive",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;
    }
    
    /* p3m parameter: mesh_off */
    else if (ARG0_IS_S("mesh_off")) {
      
      if(argc < 4) {
	Tcl_AppendResult(interp, argv[0], " needs 3 parameters",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if ((! ARG_IS_D(1, d1)) ||
	  (! ARG_IS_D(2, d2)) ||
	  (! ARG_IS_D(3, d3)))
	{
	  Tcl_AppendResult(interp, argv[0], " needs 3 DOUBLE parameters",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      if (Dp3m_set_mesh_offset(d1, d2 ,d3) == TCL_ERROR)
	{
	  Tcl_AppendResult(interp, argv[0], " parameters have to be between 0.0 an 1.0",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      argc -= 4;
      argv += 4;
    }
    
    /* p3m parameter: epsilon */
    else if(ARG0_IS_S( "epsilon")) {

      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }

      if (ARG1_IS_S("metallic")) {
	d1 = P3M_EPSILON_METALLIC;
      }
      else if (! ARG1_IS_D(d1)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter or \"metallic\"",
	                 (char *) NULL);
	return TCL_ERROR;
      }
	
      if (Dp3m_set_eps(d1) == TCL_ERROR) {
        Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
                         (char *) NULL);
        return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;	    
    }
    else {
      Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

/*****************************************************************************/


void interpolate_dipole_assignment_function()
{
  double dInterpol = 0.5 / (double)p3m.Dinter;
  int i;
  long j;

      dInterpol = 0.5 / (double)p3m.Dinter;  
    if (p3m.Dinter == 0) return;

        P3M_TRACE(fprintf(stderr,"dipolar %d - interpolating (%d) the order-%d charge assignment function\n",
		       this_node,p3m.Dinter,p3m.Dcao));

         p3m.Dinter2 = 2*p3m.Dinter + 1;

          for (i=0; i < p3m.Dcao; i++) {
             /* allocate memory for interpolation array */
             Dint_caf[i] = (double *) realloc(Dint_caf[i], sizeof(double)*(2*p3m.Dinter+1));

            /* loop over all interpolation points */
              for (j=-p3m.Dinter; j<=p3m.Dinter; j++)
                    Dint_caf[i][j+p3m.Dinter] = P3M_caf(i, j*dInterpol,p3m.Dcao);
         }
}



/*****************************************************************************/


/* assign the dipoles */
void P3M_dipole_assign()
{
  Cell *cell;
  Particle *p;
  int i,c,np,j;
  /* magnetic particle counter, dipole fraction counter */
  int cp_cnt=0;
  
  
  /* prepare local FFT mesh */
    for(i=0;i<3;i++)
      for(j=0; j<Dlm.size; j++) Drs_mesh_dip[i][j] = 0.0;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if( p[i].p.dipm != 0.0) {
	P3M_assign_dipole( p[i].r.p,p[i].p.dipm, p[i].r.dip,cp_cnt);
	cp_cnt++;
      }
    }
   } 
   DP3M_shrink_wrap_dipole_grid(cp_cnt);

}


/*****************************************************************************/



#ifdef ROTATION
/* assign the torques obtained from k-space */
static void P3M_assign_torques(double prefac, int d_rs)
{
  Cell *cell;
  Particle *p;
  int i,c,np,i0,i1,i2;
  /* particle counter, charge fraction counter */
  int cp_cnt=0, cf_cnt=0;
  /* index, index jumps for Drs_mesh array */
  int q_ind;
  int q_m_off = (Dlm.dim[2] - p3m.Dcao);
  int q_s_off = Dlm.dim[2] * (Dlm.dim[1] - p3m.Dcao);

  cp_cnt=0; cf_cnt=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i=0; i<np; i++) { 
      if( (p[i].p.dipm) != 0.0 ) {
	q_ind = Dca_fmp[cp_cnt];
	for(i0=0; i0<p3m.Dcao; i0++) {
	  for(i1=0; i1<p3m.Dcao; i1++) {
	    for(i2=0; i2<p3m.Dcao; i2++) {
/*
The following line would fill the torque with the k-space electric field
(without the self-field term) [notice the minus sign!]:		  
		    p[i].f.torque[d_rs] -= prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];;
Since the torque is the dipole moment cross-product with E, we have:	
*/
              switch (d_rs) {
		case 0:	//E_x
		  p[i].f.torque[1] -= p[i].r.dip[2]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];     
		  p[i].f.torque[2] += p[i].r.dip[1]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind]; 
		  break;
		case 1:	//E_y
		  p[i].f.torque[0] += p[i].r.dip[2]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];  
		  p[i].f.torque[2] -= p[i].r.dip[0]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];  
		  break;
		case 2:	//E_z
		  p[i].f.torque[0] -= p[i].r.dip[1]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];  
		  p[i].f.torque[1] += p[i].r.dip[0]*prefac*Dca_frac[cf_cnt]*Drs_mesh[q_ind];  
	      }
	      q_ind++; 
	      cf_cnt++;
	    }
	    q_ind += q_m_off;
	  }
	  q_ind += q_s_off;
	}
	cp_cnt++;

	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: P3M  f = (%.3e,%.3e,%.3e) in dir %d add %.5f\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],d_rs,-db_fsum));
      }
    }
  }
}
#endif
/*****************************************************************************/


/* assign the dipolar forces obtained from k-space */
static void DP3M_assign_forces_dip(double prefac, int d_rs)
{
  Cell *cell;
  Particle *p;
  int i,c,np,i0,i1,i2;
  /* particle counter, charge fraction counter */
  int cp_cnt=0, cf_cnt=0;
  /* index, index jumps for Drs_mesh array */
  int q_ind;
  int q_m_off = (Dlm.dim[2] - p3m.Dcao);
  int q_s_off = Dlm.dim[2] * (Dlm.dim[1] - p3m.Dcao);

  cp_cnt=0; cf_cnt=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i=0; i<np; i++) { 
      if( (p[i].p.dipm) != 0.0 ) {
	q_ind = Dca_fmp[cp_cnt];
	for(i0=0; i0<p3m.Dcao; i0++) {
	  for(i1=0; i1<p3m.Dcao; i1++) {
	    for(i2=0; i2<p3m.Dcao; i2++) {
	      p[i].f.f[d_rs] += prefac*Dca_frac[cf_cnt]*
	                          ( Drs_mesh_dip[0][q_ind]*p[i].r.dip[0]
		                  +Drs_mesh_dip[1][q_ind]*p[i].r.dip[1]
				  +Drs_mesh_dip[2][q_ind]*p[i].r.dip[2]);
	      q_ind++;
	      cf_cnt++;
	    }
	    q_ind += q_m_off;
	  }
	  q_ind += q_s_off;
	}
	cp_cnt++;

	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: P3M  f = (%.3e,%.3e,%.3e) in dir %d add %.5f\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],d_rs,-db_fsum));
      }
    }
  }
}

/*****************************************************************************/


double P3M_calc_kspace_forces_for_dipoles(int force_flag, int energy_flag) 
{
  int i,d,d_rs,ind,j[3];
  /**************************************************************/
   /* k space energy */
  double dipole_prefac;
  double k_space_energy_dip=0.0, node_k_space_energy_dip=0.0;
  double tmp0,tmp1;

  P3M_TRACE(fprintf(stderr,"%d: dipolar p3m_perform: \n",this_node));

  dipole_prefac = coulomb.Dprefactor / (double)(p3m.Dmesh[0]*p3m.Dmesh[1]*p3m.Dmesh[2]);
 
  if (p3m_sum_mu2 > 0) { 
    /* Gather information for FFT grid inside the nodes domain (inner local mesh) */
    /* and Perform forward 3D FFT (Charge Assignment Mesh). */
    Dgather_fft_grid(Drs_mesh_dip[0]);
    Dgather_fft_grid(Drs_mesh_dip[1]);
    Dgather_fft_grid(Drs_mesh_dip[2]);
    Dfft_perform_forw(Drs_mesh_dip[0]);
    Dfft_perform_forw(Drs_mesh_dip[1]);
    Dfft_perform_forw(Drs_mesh_dip[2]);
    //Note: after these calls, the grids are in the order yzx and not xyz anymore!!!
  }
  
  /* === K Space Calculations === */
  P3M_TRACE(fprintf(stderr,"%d: dipolar p3m_perform: k-Space\n",this_node));

  /* === K Space Energy Calculation  === */
  if(energy_flag) {
/*********************
   Dipolar energy
**********************/
  if (p3m_sum_mu2 > 0) {
    P3M_TRACE(fprintf(stderr,"%d: dipolar p3m start Energy calculation: k-Space\n",this_node));


    /* i*k differentiation for dipolar gradients: |(\Fourier{\vect{mu}}(k)\cdot \vect{k})|^2 */
    ind=0;
    i=0;
    for(j[0]=0; j[0]<Dfft_plan[3].new_mesh[0]; j[0]++) {
      for(j[1]=0; j[1]<Dfft_plan[3].new_mesh[1]; j[1]++) {
	for(j[2]=0; j[2]<Dfft_plan[3].new_mesh[2]; j[2]++) {	 
	  node_k_space_energy_dip += Dg_energy[i] * (
	  SQR(Drs_mesh_dip[0][ind]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
	      Drs_mesh_dip[1][ind]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
	      Drs_mesh_dip[2][ind]*Dd_op[j[1]+Dfft_plan[3].start[2]]
	  ) +
	  SQR(Drs_mesh_dip[0][ind+1]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
	      Drs_mesh_dip[1][ind+1]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
	      Drs_mesh_dip[2][ind+1]*Dd_op[j[1]+Dfft_plan[3].start[2]]
	      ));
	  ind += 2;
	  i++;
	}
      }
    }
    node_k_space_energy_dip *= dipole_prefac * PI / box_l[0];
    MPI_Reduce(&node_k_space_energy_dip, &k_space_energy_dip, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   
   if (Dflag_constants_energy_dipolar==0) {Dflag_constants_energy_dipolar=1;}
   compute_constants_energy_dipolar(); 
   
    k_space_energy_dip+= coulomb.Dprefactor*Dipolar_energy_correction; /* add the dipolar energy correction due to systematic Madelung-Self effects */  
   
   /*fprintf(stderr,"p3m.Depsilon=%lf\n", p3m.Depsilon);
   fprintf(stderr,"*Dipolar_energy_correction=%20.15lf\n",Dipolar_energy_correction);*/
   
    if(this_node==0) {
      /* self energy correction */
      k_space_energy_dip -= coulomb.Dprefactor*(p3m_sum_mu2*2*pow(p3m.Dalpha_L*box_l_i[0],3) * wupii/3.0);
    }

    P3M_TRACE(fprintf(stderr,"%d: dipolar p3m end Energy calculation: k-Space\n",this_node));

}
} //if (energy_flag)

  /* === K Space Force Calculation  === */
  if(force_flag) {
  /***************************        
   DIPOLAR TORQUES (k-space)
****************************/
  if (p3m_sum_mu2 > 0) {
 #ifdef ROTATION
   P3M_TRACE(fprintf(stderr,"%d: dipolar p3m start torques calculation: k-Space\n",this_node));

    /* fill in ks_mesh array for torque calculation */
    ind=0;
    i=0;
       
    for(j[0]=0; j[0]<Dfft_plan[3].new_mesh[0]; j[0]++) {	     //j[0]=n_y
      for(j[1]=0; j[1]<Dfft_plan[3].new_mesh[1]; j[1]++) {    //j[1]=n_z
	for(j[2]=0; j[2]<Dfft_plan[3].new_mesh[2]; j[2]++) {  //j[2]=n_x
	  //tmp0 = Re(mu)*k,   tmp1 = Im(mu)*k
	  
	  tmp0 = Drs_mesh_dip[0][ind]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
		 Drs_mesh_dip[1][ind]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
		 Drs_mesh_dip[2][ind]*Dd_op[j[1]+Dfft_plan[3].start[2]];
		 
	  tmp1 = Drs_mesh_dip[0][ind+1]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
		 Drs_mesh_dip[1][ind+1]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
		 Drs_mesh_dip[2][ind+1]*Dd_op[j[1]+Dfft_plan[3].start[2]];
		 
	  /* the optimal influence function is the same for torques
	     and energy */ 
	     
 	  Dks_mesh[ind]   = tmp0*Dg_energy[i]; 
	  Dks_mesh[ind+1] = tmp1*Dg_energy[i];
	  ind += 2;
	  i++;
	}
      }
    }
 
        
    /* Force component loop */
    for(d=0;d<3;d++) {
      d_rs = (d+Dks_pnum)%3;
      ind=0;
      for(j[0]=0; j[0]<Dfft_plan[3].new_mesh[0]; j[0]++) {
	for(j[1]=0; j[1]<Dfft_plan[3].new_mesh[1]; j[1]++) {
	  for(j[2]=0; j[2]<Dfft_plan[3].new_mesh[2]; j[2]++) {
	    Drs_mesh[ind] = Dd_op[ j[d]+Dfft_plan[3].start[d] ]*Dks_mesh[ind]; ind++;
	    Drs_mesh[ind] = Dd_op[ j[d]+Dfft_plan[3].start[d] ]*Dks_mesh[ind]; ind++;
	  }
	}
      }


      /* Back FFT force component mesh */
      Dfft_perform_back(Drs_mesh);
      /* redistribute force component mesh */
      Dspread_force_grid(Drs_mesh);  
      /* Assign force component from mesh to particle */
      P3M_assign_torques(dipole_prefac*(2*PI/box_l[0]), d_rs);
    }
 #endif  /*if def ROTATION */ 
    
/***************************
   DIPOLAR FORCES (k-space)
****************************/
    P3M_TRACE(fprintf(stderr,"%d: dipolar p3m start forces calculation: k-Space\n",this_node));

//Compute forces after torques because the algorithm below overwrites the grids Drs_mesh_dip !
//Note: I'll do here 9 inverse FFTs. By symmetry, we can reduce this number to 6 !
    /* fill in ks_mesh array for force calculation */
    ind=0;
    i=0;
    for(j[0]=0; j[0]<Dfft_plan[3].new_mesh[0]; j[0]++) {	     //j[0]=n_y
      for(j[1]=0; j[1]<Dfft_plan[3].new_mesh[1]; j[1]++) {    //j[1]=n_z
	for(j[2]=0; j[2]<Dfft_plan[3].new_mesh[2]; j[2]++) {  //j[2]=n_x
	  //tmp0 = Im(mu)*k,   tmp1 = -Re(mu)*k
	  tmp0 = Drs_mesh_dip[0][ind+1]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
		 Drs_mesh_dip[1][ind+1]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
		 Drs_mesh_dip[2][ind+1]*Dd_op[j[1]+Dfft_plan[3].start[2]];
	  tmp1 = Drs_mesh_dip[0][ind]*Dd_op[j[2]+Dfft_plan[3].start[0]]+
		 Drs_mesh_dip[1][ind]*Dd_op[j[0]+Dfft_plan[3].start[1]]+
		 Drs_mesh_dip[2][ind]*Dd_op[j[1]+Dfft_plan[3].start[2]];
//	  Dks_mesh[ind]   = tmp0*Dg[i];
//	  Dks_mesh[ind+1] = -tmp1*Dg[i];
          /* Next two lines modified by JJCP 26-4-06 */
	  Dks_mesh[ind]   = tmp0*Dg_force[i];
	  Dks_mesh[ind+1] = -tmp1*Dg_force[i];
	  ind += 2;
	  i++;
	}
      }
    }

    /* Force component loop */
    for(d=0;d<3;d++) {       /* direction in k space: */
    d_rs = (d+Dks_pnum)%3;
    ind=0;
    for(j[0]=0; j[0]<Dfft_plan[3].new_mesh[0]; j[0]++) {	     //j[0]=n_y
      for(j[1]=0; j[1]<Dfft_plan[3].new_mesh[1]; j[1]++) {    //j[1]=n_z
	for(j[2]=0; j[2]<Dfft_plan[3].new_mesh[2]; j[2]++) {  //j[2]=n_x
	  tmp0 = Dd_op[ j[d]+Dfft_plan[3].start[d] ]*Dks_mesh[ind];
	  Drs_mesh_dip[0][ind] = Dd_op[ j[2]+Dfft_plan[3].start[d] ]*tmp0;
	  Drs_mesh_dip[1][ind] = Dd_op[ j[0]+Dfft_plan[3].start[d] ]*tmp0;
	  Drs_mesh_dip[2][ind] = Dd_op[ j[1]+Dfft_plan[3].start[d] ]*tmp0;
	  ind++;
	  tmp0 = Dd_op[ j[d]+Dfft_plan[3].start[d] ]*Dks_mesh[ind];
	  Drs_mesh_dip[0][ind] = Dd_op[ j[2]+Dfft_plan[3].start[d] ]*tmp0;
	  Drs_mesh_dip[1][ind] = Dd_op[ j[0]+Dfft_plan[3].start[d] ]*tmp0;
	  Drs_mesh_dip[2][ind] = Dd_op[ j[1]+Dfft_plan[3].start[d] ]*tmp0;
	  ind++;
	}
      }
    }

      /* Back FFT force component mesh */
      Dfft_perform_back(Drs_mesh_dip[0]);
      Dfft_perform_back(Drs_mesh_dip[1]);
      Dfft_perform_back(Drs_mesh_dip[2]);
      /* redistribute force component mesh */
      Dspread_force_grid(Drs_mesh_dip[0]);
      Dspread_force_grid(Drs_mesh_dip[1]);
      Dspread_force_grid(Drs_mesh_dip[2]);
      /* Assign force component from mesh to particle */
      DP3M_assign_forces_dip(dipole_prefac*pow(2*PI/box_l[0],2), d_rs); 
   }
   
       P3M_TRACE(fprintf(stderr,"%d: dipolar p3m end forces calculation: k-Space\n",this_node));

   
 } /* of if (p3m_sum_mu2>0 */
} /* of if(force_flag) */

 
  if (p3m.Depsilon != P3M_EPSILON_METALLIC) {
    k_space_energy_dip += calc_surface_term(force_flag, energy_flag);
   }


  return k_space_energy_dip;
}



/************************************************************/

double calc_surface_term(int force_flag, int energy_flag)
{
 
   int np, c, i,ip;
  Particle *part;
  double pref =coulomb.Dprefactor*4*M_PI*box_l_i[0]*box_l_i[1]*box_l_i[2]/(2*p3m.Depsilon + 1);
  double suma,ax,ay,az;
  double en;
  double  *sumix=NULL,*sumiy=NULL,*sumiz=NULL;
  double  *mx=NULL,*my=NULL,*mz=NULL;
     
      
     
  if(n_nodes==1) {
 
 
     // We put all the dipolar momenta in a the arrays mx,my,mz according to the id-number of the particles   
     mx = (double *) malloc(sizeof(double)*n_total_particles);
     my = (double *) malloc(sizeof(double)*n_total_particles);
     mz = (double *) malloc(sizeof(double)*n_total_particles);
    
  
     for (c = 0; c < local_cells.n; c++) {
       np   = local_cells.cell[c]->n;
       part = local_cells.cell[c]->part;
       for (i = 0; i < np; i++){
 	 ip=part[i].p.identity;
	 mx[ip]=part[i].r.dip[0];
	 my[ip]=part[i].r.dip[1];
	 mz[ip]=part[i].r.dip[2];	 
      }  
     } 



     // we will need the sum of all dipolar momenta vectors    
      ax=0.0;
      ay=0.0;
      az=0.0;

      for (i = 0; i < n_total_particles; i++){
         ax+=mx[i];
         ay+=my[i];
         az+=mz[i];
      }   

 
     
     //for (i = 0; i < n_total_particles; i++){
     //  fprintf(stderr,"part ip:%d, mux: %le, muy:%le, muz:%le \n",i,mx[i],my[i],mz[i]);   
     //}
     
 
     //Now we can proceed to compute things .....
  
     if (energy_flag) {
      
        suma=0.0;
        for (i = 0; i < n_total_particles; i++){
 	      suma+=mx[i]*ax+my[i]*ay+mz[i]*az;
        }  	      

        //fprintf(stderr,"energia, pref=%le, suma=%le\n",pref,suma);

        en = 0.5*pref*suma;
       
     } else {
        en = 0;
     } 
      

      
       
     if (force_flag) {
 
 
          //fprintf(stderr," number of particles= %d ",n_total_particles);   

          sumix = (double *) malloc(sizeof(double)*n_total_particles);
          sumiy = (double *) malloc(sizeof(double)*n_total_particles);
          sumiz = (double *) malloc(sizeof(double)*n_total_particles);
	  
          for (i = 0; i < n_total_particles; i++){
	    sumix[i]=my[i]*az-mz[i]*ay;
            sumiy[i]=mz[i]*ax-mx[i]*az;
            sumiz[i]=mx[i]*ay-my[i]*ax;
	  }
	  
   
   
         // for (i = 0; i < n_total_particles; i++){
  	 //    fprintf(stderr,"part %d, correccions torque  x:%le, y:%le, z:%le\n",i,sumix[i],sumiy[i],sumiz[i]);
         // }
	      
          #ifdef ROTATION	      
          for (c = 0; c < local_cells.n; c++) {
             np	= local_cells.cell[c]->n;
             part = local_cells.cell[c]->part;
             for (i = 0; i < np; i++){
	     
	        ip=part[i].p.identity;
		
	        // fprintf(stderr,"part %d, torque abans %le %le %le\n",ip,part[i].f.torque[0],part[i].f.torque[1],part[i].f.torque[2]);
	      
		part[i].f.torque[0] -= pref*sumix[ip];
		part[i].f.torque[1] -= pref*sumiy[ip];
		part[i].f.torque[2] -= pref*sumiz[ip];
		
	       // fprintf(stderr,"part %d, torque despres %le %le %le\n",ip,part[i].f.torque[0],part[i].f.torque[1],part[i].f.torque[2]);
 	     }	
          }
          #endif
	     
	  free(sumix);     
  	  free(sumiy);     
	  free(sumiz);     
     }
       
    free(mx);	 
    free(my);	 
    free(mz);	 
    	
  } else {
        //The code is not prepared to run in more than one node
      
      fprintf(stderr,"dipolar-P3M: Non metallic environment is  not ready to work in more than one Node, Sorry ....");
      fprintf(stderr," I am NOT computing the surface term for the dipolar interaction.neither fortorques not for the energy...");
      return 0.;

  }
 
 	    
  return en;
 
}


/************************************************************/
void Dgather_fft_grid(double* themesh)
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;

  P3M_TRACE(fprintf(stderr,"%d: Dgather_fft_grid:\n",this_node));

  /* direction loop */
  for(s_dir=0; s_dir<6; s_dir++) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(Dsm.s_size[s_dir]>0) 
      pack_block(themesh, Dsend_grid, Dsm.s_ld[s_dir], Dsm.s_dim[s_dir], Dlm.dim, 1);
      
    /* communication */
    if(node_neighbors[s_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[s_dir/2]+evenodd)%2==0) {
	  if(Dsm.s_size[s_dir]>0) 
	    MPI_Send(Dsend_grid, Dsm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], REQ_P3M_GATHER_D, MPI_COMM_WORLD);
	}
	else {
	  if(Dsm.r_size[r_dir]>0) 
	    MPI_Recv(Drecv_grid, Dsm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], REQ_P3M_GATHER_D, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = Drecv_grid;
      Drecv_grid = Dsend_grid;
      Dsend_grid = tmp_ptr;
    }
    /* add recv block */
    if(Dsm.r_size[r_dir]>0) {
      add_block(Drecv_grid, themesh, Dsm.r_ld[r_dir], Dsm.r_dim[r_dir], Dlm.dim); 
    }
  }
}



/************************************************************/


void Dspread_force_grid(double* themesh)
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;
  P3M_TRACE(fprintf(stderr,"%d: dipolar spread_force_grid:\n",this_node));

  /* direction loop */
  for(s_dir=5; s_dir>=0; s_dir--) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(Dsm.s_size[s_dir]>0) 
      pack_block(themesh, Dsend_grid, Dsm.r_ld[r_dir], Dsm.r_dim[r_dir], Dlm.dim, 1);
    /* communication */
    if(node_neighbors[r_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[r_dir/2]+evenodd)%2==0) {
	  if(Dsm.r_size[r_dir]>0) 
	    MPI_Send(Dsend_grid, Dsm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], REQ_P3M_SPREAD_D, MPI_COMM_WORLD);
   	}
	else {
	  if(Dsm.s_size[s_dir]>0) 
	    MPI_Recv(Drecv_grid, Dsm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], REQ_P3M_SPREAD_D, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = Drecv_grid;
      Drecv_grid = Dsend_grid;
      Dsend_grid = tmp_ptr;
    }
    /* un pack recv block */
    if(Dsm.s_size[s_dir]>0) {
      unpack_block(Drecv_grid, themesh, Dsm.s_ld[s_dir], Dsm.s_dim[s_dir], Dlm.dim, 1); 
    }
  }
}


/*****************************************************************************/

void Drealloc_ca_fields(int newsize)
{
  newsize = ((newsize + CA_INCREMENT - 1)/CA_INCREMENT)*CA_INCREMENT;
  if (newsize == Dca_num) return;
  if (newsize < CA_INCREMENT) newsize = CA_INCREMENT;

   P3M_TRACE(fprintf(stderr,"%d: realloc_ca_fields: dipolar,  old_size=%d -> new_size=%d\n",this_node,Dca_num,newsize));
   Dca_num = newsize;
   Dca_frac = (double *)realloc(Dca_frac, p3m.Dcao3*Dca_num*sizeof(double));
   Dca_fmp  = (int *)realloc(Dca_fmp, Dca_num*sizeof(int));
  
}


/*****************************************************************************/


void Dcalc_meshift(void)
{
  int i;
  double dmesh;
     dmesh = (double)p3m.Dmesh[0];
     Dmeshift = (double *) realloc(Dmeshift, p3m.Dmesh[0]*sizeof(double));
     for (i=0; i<p3m.Dmesh[0]; i++) Dmeshift[i] = i - dround(i/dmesh)*dmesh;   
}



/*****************************************************************************/


void Dcalc_differential_operator()
{
  int i;
  double dmesh;

  dmesh = (double)p3m.Dmesh[0];
  Dd_op = (double *) realloc(Dd_op, p3m.Dmesh[0]*sizeof(double));

  for (i=0; i<p3m.Dmesh[0]; i++) 
    Dd_op[i] = (double)i - dround((double)i/dmesh)*dmesh;

    Dd_op[p3m.Dmesh[0]/2] = 0;
}

/*****************************************************************************/


void Dcalc_influence_function_force()
{
  int i,n[3],ind;
  int end[3];
  int size=1;
  double fak1,fak2;
  double nominator[1]={0.0},denominator=0.0;

  Dcalc_meshift();

  for(i=0;i<3;i++) {
    size *= Dfft_plan[3].new_mesh[i];
    end[i] = Dfft_plan[3].start[i] + Dfft_plan[3].new_mesh[i];
  }
  Dg_force = (double *) realloc(Dg_force, size*sizeof(double));
  fak1  = p3m.Dmesh[0]*p3m.Dmesh[0]*p3m.Dmesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=Dfft_plan[3].start[0]; n[0]<end[0]; n[0]++)
    for(n[1]=Dfft_plan[3].start[1]; n[1]<end[1]; n[1]++)
      for(n[2]=Dfft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-Dfft_plan[3].start[2]) + Dfft_plan[3].new_mesh[2] * ((n[1]-Dfft_plan[3].start[1]) + (Dfft_plan[3].new_mesh[1]*(n[0]-Dfft_plan[3].start[0])));

	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  Dg_force[ind] = 0.0;
	else if( (n[0]%(p3m.Dmesh[0]/2)==0) &&
		 (n[1]%(p3m.Dmesh[0]/2)==0) &&
		 (n[2]%(p3m.Dmesh[0]/2)==0) )
	  Dg_force[ind] = 0.0;
	else {
	  denominator = Dperform_aliasing_sums_force(n,nominator);
	  fak2 =  nominator[0];
	  fak2 /= pow(SQR(Dd_op[n[0]])+SQR(Dd_op[n[1]])+SQR(Dd_op[n[2]]),3)  * SQR(denominator) ;
	  Dg_force[ind] = fak1*fak2;
	}
      }
}


/*****************************************************************************/

MDINLINE double Dperform_aliasing_sums_force(int n[3], double nominator[1])
{
  double denominator=0.0;
  /* lots of temporary variables... */
  double sx,sy,sz,f1,f2,f3,mx,my,mz,nmx,nmy,nmz,nm2,expo;
  double limit = 30;
  double n_nm;
  double n_nm3;

  nominator[0]=0.0;
  
  f1 = 1.0/(double)p3m.Dmesh[0];
  f2 = SQR(PI/(p3m.Dalpha_L));

  for(mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = Dmeshift[n[0]] + p3m.Dmesh[0]*mx;
    sx  = pow(sinc(f1*nmx),2.0*p3m.Dcao);
    for(my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = Dmeshift[n[1]] + p3m.Dmesh[0]*my;
      sy  = sx*pow(sinc(f1*nmy),2.0*p3m.Dcao);
      for(mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
	nmz = Dmeshift[n[2]] + p3m.Dmesh[0]*mz;
	sz  = sy*pow(sinc(f1*nmz),2.0*p3m.Dcao);
	
	nm2          =  SQR(nmx)+SQR(nmy)+SQR(nmz);
	expo         =  f2*nm2;
	f3           =  (expo<limit) ? sz*exp(-expo)/nm2 : 0.0;

	n_nm = Dd_op[n[0]]*nmx + Dd_op[n[1]]*nmy + Dd_op[n[2]]*nmz;
	n_nm3 = n_nm*n_nm*n_nm; 
	
	nominator[0] += f3*n_nm3;
	denominator  += sz;
      }
    }
  }
  return denominator;
}


/*****************************************************************************/

void Dcalc_influence_function_energy()
{
  int i,n[3],ind;
  int end[3];
  int size=1;
  double fak1,fak2;
  double nominator[1]={0.0},denominator=0.0;

  Dcalc_meshift();

  for(i=0;i<3;i++) {
    size *= Dfft_plan[3].new_mesh[i];
    end[i] = Dfft_plan[3].start[i] + Dfft_plan[3].new_mesh[i];
  }
  Dg_energy = (double *) realloc(Dg_energy, size*sizeof(double));
  fak1  = p3m.Dmesh[0]*p3m.Dmesh[0]*p3m.Dmesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=Dfft_plan[3].start[0]; n[0]<end[0]; n[0]++)
    for(n[1]=Dfft_plan[3].start[1]; n[1]<end[1]; n[1]++)
      for(n[2]=Dfft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-Dfft_plan[3].start[2]) + Dfft_plan[3].new_mesh[2] * ((n[1]-Dfft_plan[3].start[1]) + (Dfft_plan[3].new_mesh[1]*(n[0]-Dfft_plan[3].start[0])));

	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  Dg_energy[ind] = 0.0;
	else if( (n[0]%(p3m.Dmesh[0]/2)==0) &&
		 (n[1]%(p3m.Dmesh[0]/2)==0) &&
		 (n[2]%(p3m.Dmesh[0]/2)==0) )
	  Dg_energy[ind] = 0.0;
	else {
	  denominator = Dperform_aliasing_sums_energy(n,nominator);
	  fak2 =  nominator[0];
	  fak2 /= pow(SQR(Dd_op[n[0]])+SQR(Dd_op[n[1]])+SQR(Dd_op[n[2]]),2)  * SQR(denominator) ;
	  Dg_energy[ind] = fak1*fak2;
	}
      }
}

/*****************************************************************************/

MDINLINE double Dperform_aliasing_sums_energy(int n[3], double nominator[1])
{ 
  double denominator=0.0;
  /* lots of temporary variables... */
  double sx,sy,sz,f1,f2,f3,mx,my,mz,nmx,nmy,nmz,nm2,expo;
  double limit = 30;
  double n_nm;
  double n_nm2;

  nominator[0]=0.0;
    
  f1 = 1.0/(double)p3m.Dmesh[0];
  f2 = SQR(PI/(p3m.Dalpha_L));

  for(mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = Dmeshift[n[0]] + p3m.Dmesh[0]*mx;
    sx  = pow(sinc(f1*nmx),2.0*p3m.Dcao);
    for(my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = Dmeshift[n[1]] + p3m.Dmesh[0]*my;
      sy  = sx*pow(sinc(f1*nmy),2.0*p3m.Dcao);
      for(mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
	nmz = Dmeshift[n[2]] + p3m.Dmesh[0]*mz;
	sz  = sy*pow(sinc(f1*nmz),2.0*p3m.Dcao);
	
	nm2          =  SQR(nmx)+SQR(nmy)+SQR(nmz);
	expo         =  f2*nm2;
	f3           =  (expo<limit) ? sz*exp(-expo)/nm2 : 0.0;

	n_nm = Dd_op[n[0]]*nmx + Dd_op[n[1]]*nmy + Dd_op[n[2]]*nmz;  /* JJCP 26-4-06 */
	n_nm2 = n_nm*n_nm; 
	nominator[0] += f3*n_nm2;
	denominator  += sz;
      }
    }
  }
  return denominator;
}


/*****************************************************************************/


/************************************************
 * Functions for dipoloar P3M Parameter tuning
 * This tuning is based on the P3M tunning of the charges
 which in turn is based on the P3M_tune by M. Deserno
 ************************************************/

#define P3M_TUNE_MAX_CUTS 50

int DP3M_tune_parameters(Tcl_Interp *interp)
{
  int i,ind, try=0, best_try=0, n_cuts;
  double r_cut_iL, r_cut_iL_min  , r_cut_iL_max, r_cut_iL_best=0, cuts[P3M_TUNE_MAX_CUTS], cut_start;
  int    mesh    , mesh_min      , mesh_max    , mesh_best=0;
  int    cao     , cao_min       , cao_max     , cao_best=0;
  double alpha_L , alpha_L_best=0, accuracy    , accuracy_best=0;
  double mesh_size, k_cut;
  double rs_err, rs_err_best=0, ks_err, ks_err_best=0;
  double int_time=0, min_time=1e20, int_num;
  char b1[TCL_DOUBLE_SPACE + 12],b2[TCL_DOUBLE_SPACE + 12],b3[TCL_DOUBLE_SPACE + 12];
 
  P3M_TRACE(fprintf(stderr,"%d: DP3M_tune_parameters\n",this_node));
  
  /* preparation */
  mpi_bcast_event(P3M_COUNT_DIPOLES);

  /* calculate r_cut_iL tune range */
  if(p3m.Dr_cut_iL == 0.0) { 
    n_cuts = P3M_TUNE_MAX_CUTS;
    for(i=0;i<n_cuts;i++) {
      if(min_local_box_l == min_box_l)
	cuts[i] = min_local_box_l/(i+2.0)-(skin);
      else 
	cuts[i] = min_local_box_l/(i+1.0)-(skin);
      cuts[i]*=box_l_i[0];
      if( cuts[i] <= 0.0 ) {
	n_cuts = i; 
	break;
      } 
    }
    r_cut_iL_max = cuts[0];
    r_cut_iL_min = cuts[n_cuts-1];
  }
  else { 
    n_cuts = 1;
    r_cut_iL_min = r_cut_iL_max = p3m.Dr_cut_iL; 
    cuts[0] = p3m.Dr_cut_iL;
  }
  /* calculate mesh tune range */
  if(p3m.Dmesh[0] == 0 ) {
    double expo;
     expo = log(pow((double)p3m_sum_dip_part,(1.0/3.0)))/log(2.0);    
    mesh_min = (int)(pow(2.0,(double)((int)expo))+0.1);
    mesh_max = mesh_min*4;
    if(mesh_min < 8) { mesh_min = 8; mesh_max = 16; }
  }
  else { mesh_min = mesh_max = p3m.Dmesh[0]; }
  /* calculate cao tune range */
  if(p3m.Dcao == 0) { cao_min = 1; cao_max = 7; }
  else             { cao_min = cao_max = p3m.Dcao; }

  /* Print Status */
  sprintf(b1,"%.5e",p3m.Daccuracy);
  Tcl_AppendResult(interp, "P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);

      sprintf(b2,"%d",p3m_sum_dip_part);   

  Tcl_PrintDouble(interp, p3m_sum_mu2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, r_cut_iL_min, b1);  Tcl_PrintDouble(interp, r_cut_iL_max, b2);
  Tcl_AppendResult(interp, "Range for p3m.Dr_cut_iL: [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",mesh_min);  sprintf(b2,"%d",mesh_max);
  Tcl_AppendResult(interp, "Range for p3m.Dmesh:     [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",cao_min);  sprintf(b2,"%d",cao_max);
  Tcl_AppendResult(interp, "Range for p3m.Dcao:      [",b1,"-",b2,"]","\n\n", (char *) NULL);
  Tcl_AppendResult(interp, "set mesh cao r_cut_iL     alpha_L      err          ks_err     rs_err     time [ms]\n", (char *) NULL);

  /* Tuning Loops */
  for(mesh = mesh_min; mesh <= mesh_max; mesh*=2) { /* mesh loop */
    cut_start = box_l[0] * box_l_i[0];

       if(mesh <= 32 || p3m_sum_dip_part > 2000) int_num=5; else int_num=1;  

    for(cao = cao_min; cao <= cao_max; cao++) {     /* cao loop */
      mesh_size = box_l[0]/(double)mesh;
      k_cut =  mesh_size*cao/2.0;
      if(cao < mesh && k_cut < dmin(min_box_l,min_local_box_l)-skin) {
	ind=0;
	for(i=0;i< n_cuts -1;i++) {
	  if(cut_start <= cuts[i]) ind=i+1;
	}
	while (ind < n_cuts) {           /* r_cut_iL loop */
	  r_cut_iL = cuts[ind];
	  /* calc maximal real space error for setting */
	  //Beginning of JJCP modification 15/5/2006 
	  
	     if(r_cut_iL *box_l[0] < 1.0) break;   // It has little sense checking so little rcuts ...  
	  
	         
		 
	     //Alpha cannot be zero in the dipolar case because real_space formula breaks down	 
	     
	     // Here we follow a method different from Coulombic case because the formula for the real space error
	     // is a trascendental equation for alpha in difference to the coulomb case 
	     
	    
	     
	    rs_err=P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,0.001);
	     

	  if(sqrt(2.0)*rs_err > p3m.Daccuracy) {
	    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
	    alpha_L = sqrt(log(sqrt(2.0)*rs_err/p3m.Daccuracy)) / r_cut_iL;
	    /* calculate real space and k space error for this alpha_L */

	    alpha_L=JJ_rtbisection(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,
                    0.0001*box_l[0],5.0*box_l[0],0.0001,p3m.Daccuracy);
	    
	    
	    rs_err = P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,alpha_L);
	    ks_err = P3M_DIPOLAR_k_space_error(box_l[0],coulomb.Dprefactor,mesh,cao,p3m_sum_dip_part,p3m_sum_mu2,alpha_L);

	     accuracy = sqrt(SQR(rs_err)+SQR(ks_err));
	    
	    /* check if this matches the accuracy goal */
	    if(accuracy <= p3m.Daccuracy) {
	      cut_start = cuts[ind];
	      /* broadcast p3m parameters for test run */
	      p3m.Dr_cut_iL = r_cut_iL;
	      p3m.Dmesh[0]  = p3m.Dmesh[1] = p3m.Dmesh[2] = mesh;
	      p3m.Dcao      = cao;
	      p3m.Dalpha_L  = alpha_L;
	      P3M_scaleby_box_l_dipoles();
	      /* initialize p3m structures */
	      mpi_bcast_coulomb_params();
	      /* perform force calculation test */
	      int_time = time_force_calc(int_num);
	      if (int_time == -1)
		return TCL_ERROR;
	      try++;
	      P3M_TRACE(fprintf(stderr,"%d ",try));
	      /* print result */
	      sprintf(b1,"%-3d",try); sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
	      Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
	      sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",alpha_L); sprintf(b3,"%.5e",accuracy);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
	      sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err); sprintf(b3,"%-8d",(int)int_time);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
	      if(int_time <= min_time  && r_cut_iL > 0) {
		min_time      = int_time;
		r_cut_iL_best = r_cut_iL;
		mesh_best     = mesh;
		cao_best      = cao;
		alpha_L_best  = alpha_L;
		accuracy_best = sqrt(SQR(rs_err)+SQR(ks_err));
		rs_err_best   = rs_err;
		ks_err_best   = ks_err;
		best_try      = try;
	      }
	    }
	  }
	  ind++;
	}
      }
    }
  }
  P3M_TRACE(fprintf(stderr,"\n"));
  if(try==0) {
    Tcl_AppendResult(interp, "\nFailed to tune P3M parameters to required accuracy ! \n", (char *) NULL);
    return (TCL_ERROR);
  }

  /* set tuned p3m parameters */
  p3m.Dr_cut_iL = r_cut_iL_best;
  p3m.Dmesh[0]  = p3m.Dmesh[1] = p3m.Dmesh[2] = mesh_best;
  p3m.Dcao      = cao_best;
  p3m.Dalpha_L  = alpha_L_best;
  p3m.Daccuracy = accuracy_best;
  P3M_scaleby_box_l_dipoles();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  sprintf(b1,"%d",try);
  Tcl_AppendResult(interp, "\nTune results of ",b1," trials:\n", (char *) NULL);
  sprintf(b1,"%-3d",best_try); sprintf(b2,"%-4d",mesh_best); sprintf(b3,"%-3d",cao_best);
  Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL_best); sprintf(b2,"%.5e",alpha_L_best); sprintf(b3,"%.5e",accuracy_best);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
  sprintf(b1,"%.3e",rs_err_best); sprintf(b2,"%.3e",ks_err_best); sprintf(b3,"%-8d",(int)min_time);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
  sprintf(b1,"%g",coulomb.Dbjerrum); sprintf(b2,"%g",p3m.Dr_cut); sprintf(b3,"%d",mesh_best); 
  Tcl_AppendResult(interp, "=> inter coulomb ", b1, " p3m ", b2, " ", b3, (char *) NULL);
  sprintf(b1,"%d",cao_best); sprintf(b2,"%g",p3m.Dalpha); sprintf(b3,"%g",accuracy_best);
  Tcl_AppendResult(interp, " ", b1," ", b2," ", b3," \n", (char *) NULL);

  return (TCL_OK);
}

/*****************************************************************************/


/** get the minimal error for this combination of parameters. In fact, the real space error is tuned such that it
    contributes half of the total error, and then the Fourier space error is calculated. Returns the error and the
    optimal alpha, or 0 if this combination does not work at all */
static double Dget_accuracy(int mesh, int cao, double r_cut_iL, double *_alpha_L, double *_rs_err, double *_ks_err)
{
  double rs_err, ks_err;
  double alpha_L;
  P3M_TRACE(fprintf(stderr, "Dget_accuracy: mesh %d, cao %d, r_cut %f ", mesh, cao, r_cut_iL));

  /* calc maximal real space error for setting */

    //Alpha cannot be zero in the dipolar case because real_space formula breaks down	     
    //Idem of the previous function DP3M_tune_parameters, here we do nothing
    rs_err =P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,0.001);
    
  
    if(M_SQRT2*rs_err > p3m.Daccuracy) {
     /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
         alpha_L=JJ_rtbisection(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,
   	     0.0001*box_l[0],5.0*box_l[0],0.0001,p3m.Daccuracy);

    }

  else
    /* even alpha=0 is ok, however, we cannot choose it since it kills the k-space error formula.
       Anyways, this very likely NOT the optimal solution */
    alpha_L = 0.1;

  *_alpha_L = alpha_L;
  /* calculate real space and k space error for this alpha_L */

    rs_err = P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,p3m_sum_dip_part,p3m_sum_mu2,alpha_L);
    ks_err = P3M_DIPOLAR_k_space_error(box_l[0],coulomb.Dprefactor,mesh,cao,p3m_sum_dip_part,p3m_sum_mu2,alpha_L);

  *_rs_err = rs_err;
  *_ks_err = ks_err;
  P3M_TRACE(fprintf(stderr, "dipolar tuning resulting: %f -> %f %f\n", alpha_L, rs_err, ks_err));
  return sqrt(SQR(rs_err)+SQR(ks_err));
}


/*****************************************************************************/

/** get the optimal alpha and the corresponding computation time for fixed mesh, cao, r_cut and alpha */
static double Dp3m_mcr_time(int mesh, int cao, double r_cut_iL, double alpha_L)
{
  /* rounded up 2000/n_charges timing force evaluations */
    int int_num = (1999 + p3m_sum_dip_part)/p3m_sum_dip_part;

  /* broadcast p3m parameters for test run */
  p3m.Dr_cut_iL = r_cut_iL;
  p3m.Dmesh[0]  = p3m.Dmesh[1] = p3m.Dmesh[2] = mesh;
  p3m.Dcao      = cao;
  p3m.Dalpha_L  = alpha_L;
  P3M_scaleby_box_l_dipoles();
  /* initialize p3m structures */
  mpi_bcast_coulomb_params();
  /* perform force calculation test */
  return time_force_calc(int_num);    
}


/*****************************************************************************/
/** get the optimal alpha and the corresponding computation time for fixed mesh, cao. The r_cut is determined via
    a simple bisection. Returns -1 if the force evaluation does not work, -2 if there is no valid r_cut, and -3 if
    the charge assigment order is to large for this grid */
static double Dp3m_mc_time(Tcl_Interp *interp, int mesh, int cao,
			  double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			  double *_alpha_L, double *_accuracy)
{
  double int_time;
  double r_cut_iL;
  double rs_err, ks_err, mesh_size, k_cut;
  int i, n_cells;
  char b1[TCL_DOUBLE_SPACE + 12],b2[TCL_DOUBLE_SPACE + 12],b3[TCL_DOUBLE_SPACE + 12];
  /* initial checks. */
  mesh_size = box_l[0]/(double)mesh;
  k_cut =  mesh_size*cao/2.0;
  P3M_TRACE(fprintf(stderr, "Dp3m_mc_time: mesh=%d, cao=%d, rmin=%f, rmax=%f\n",
		    mesh, cao, r_cut_iL_min, r_cut_iL_max));
  if(cao >= mesh || k_cut >= dmin(min_box_l,min_local_box_l) - skin) {
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," cao too large for this mesh\n", (char *) NULL);
    return -3;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary fails, there is no possible r_cut */
  if ((*_accuracy = Dget_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err, &ks_err)) > p3m.Daccuracy) {
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  accuracy not achieved\n", (char *) NULL);
    return -2;
  }

  for (;;) {
    P3M_TRACE(fprintf(stderr, "Dp3m_mc_time: interval [%f,%f]\n", r_cut_iL_min, r_cut_iL_max));
    r_cut_iL = 0.5*(r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    if (Dget_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err) > p3m.Daccuracy)
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }
  /* final result is always the upper interval boundary, since only there
     we know that the desired minimal accuracy is obtained */
  *_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* check whether we are running P3M+DLC, and whether we leave a reasonable gap space */
  if (coulomb.Dmethod == DIPOLAR_MDLC_P3M)   fprintf(stderr, "tunning when dlc needs to be fixed, p3m-dipoles.c \n");
  
  /*
  needs to be fixed
  if (coulomb.method == DIPOLAR_MDLC_P3M && elc_params.gap_size <= 1.1*r_cut_iL*box_l[0]) {
    P3M_TRACE(fprintf(stderr, "Dp3m_mc_time: mesh %d cao %d r_cut %f reject r_cut %f > gap %f\n", mesh, cao, r_cut_iL,
		      2*r_cut_iL*box_l[0], elc_params.gap_size));
    // print result 
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  conflict with ELC\n", (char *) NULL);
    return -2;
  }
  */

  /* check whether this radius is too large, so that we would use less cells than allowed */
  n_cells = 1;
  for (i = 0; i < 3; i++)
    n_cells *= (int)(floor(local_box_l[i]/(r_cut_iL*box_l[0] + skin)));
  if (n_cells < min_num_cells) {
    P3M_TRACE(fprintf(stderr, "Dp3m_mc_time: mesh %d cao %d r_cut %f reject n_cells %d\n", mesh, cao, r_cut_iL, n_cells));
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  radius dangerously high\n", (char *) NULL);
    return -2;
  }

  int_time = Dp3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -1) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "tuning failed, test integration not possible", (char *)NULL);
    return -1;
  }

  *_accuracy = Dget_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);

  P3M_TRACE(fprintf(stderr, "Dp3m_mc_time: mesh %d cao %d r_cut %f time %f\n", mesh, cao, r_cut_iL, int_time));
  /* print result */
  sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
  Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
  sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err); sprintf(b3,"%-8d",(int)int_time);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);

  return int_time;
}


/*****************************************************************************/

/** get the optimal alpha and the corresponding computation time for fixed mesh. *cao
    should contain an initial guess, which is then adapted by stepping up and down. Returns the time
    upon completion, -1 if the force evaluation does not work, and -2 if the accuracy cannot be met */
static double Dp3m_m_time(Tcl_Interp *interp, int mesh,
			 int cao_min, int cao_max, int *_cao,
			 double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			 double *_alpha_L, double *_accuracy)
{
  double best_time = -1, tmp_time, tmp_r_cut_iL, tmp_alpha_L=0.0, tmp_accuracy=0.0;
  /* in which direction improvement is possible. Initially, we dont know it yet. */
  int final_dir = 0;
  int cao = *_cao;

  P3M_TRACE(fprintf(stderr, "Dp3m_m_time: Dmesh=%d, Dcao_min=%d, Dcao_max=%d, Drmin=%f, Drmax=%f\n",
		    mesh, cao_min, cao_max, r_cut_iL_min, r_cut_iL_max));
  /* the initial step sets a timing mark. If there is no valid r_cut, we can only try
     to increase cao to increase the obtainable precision of the far formula. */
  do {
    tmp_time = Dp3m_mc_time(interp, mesh, cao,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -1) return -1;
    /* cao is too large for this grid, but still the accuracy cannot be achieved, give up */
    if (tmp_time == -3) {
      P3M_TRACE(fprintf(stderr, "Dp3m_m_time: no possible cao found\n"));
      return -2;
    }
    /* we have a valid time, start optimising from there */
    if (tmp_time >= 0) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
      break;
    }
    /* the required accuracy could not be obtained, try higher caos. Therefore optimisation can only be
       obtained with even higher caos, but not lower ones */
    P3M_TRACE(fprintf(stderr, "Dp3m_m_time: doesn't give precision, step up\n"));
    cao++;
    final_dir = 1;
  }
  while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max) return -2;

  /* at the boundaries, only the opposite direction can be used for optimisation */
  if (cao == cao_min)      final_dir = 1;
  else if (cao == cao_max) final_dir = -1;

  P3M_TRACE(fprintf(stderr, "Dp3m_m_time: final constraints dir %d\n", final_dir));

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
	Dp3m_mc_time(interp, mesh, cao + final_dir,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
      /* bail out on errors, as usual */
      if (tmp_time == -1) return -1;
      /* in this direction, we cannot optimise, since we get into precision trouble */
      if (tmp_time < 0) continue;

      if (tmp_time < best_time) {
	best_time  = tmp_time;
	*_r_cut_iL = tmp_r_cut_iL;
	*_alpha_L  = tmp_alpha_L;
	*_accuracy = tmp_accuracy;
	*_cao      = cao + final_dir;
      }
    }
    /* choose the direction which was optimal, if any of the two */
    if      (dir_times[0] == best_time) { final_dir = -1; }
    else if (dir_times[2] == best_time) { final_dir = 1; }
    else {
      /* no improvement in either direction, however if one is only marginally worse, we can still try*/
      /* down is possible and not much worse, while up is either illegal or even worse */
      if ((dir_times[0] >= 0 && dir_times[0] < best_time + P3M_TIME_GRAN) &&
	  (dir_times[2] < 0 || dir_times[2] > dir_times[0]))
	final_dir = -1;
      /* same for up */
      else if ((dir_times[2] >= 0 && dir_times[2] < best_time + P3M_TIME_GRAN) &&
	       (dir_times[0] < 0 || dir_times[0] > dir_times[2]))
	final_dir = 1;
      else {
	/* really no chance for optimisation */
	P3M_TRACE(fprintf(stderr, "Dp3m_m_time: Dmesh=%d final Dcao=%d time=%f\n",mesh, cao, best_time));
	return best_time;
      }
    }
    /* we already checked the initial cao and its neighbor */
    cao += 2*final_dir;
  }
  else {
    /* here some constraint is active, and we only checked the initial cao itself */
    cao += final_dir;
  }

  P3M_TRACE(fprintf(stderr, "Dp3m_m_time: optimise in direction %d\n", final_dir));

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time = Dp3m_mc_time(interp, mesh, cao,  r_cut_iL_min, r_cut_iL_max,
			   &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out on errors, as usual */
    if (tmp_time == -1) return -1;
    /* if we cannot meet the precision anymore, give up */
    if (tmp_time < 0) break;

    if (tmp_time < best_time) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
    }
    /* no hope of further optimisation */
    else if (tmp_time > best_time + P3M_TIME_GRAN)
      break;
  }
  P3M_TRACE(fprintf(stderr, "Dp3m_m_time: Dmesh=%d final Dcao=%d Dr_cut=%f time=%f\n",mesh, *_cao, *_r_cut_iL, best_time));
  return best_time;
}


/*****************************************************************************/

int DP3M_adaptive_tune_parameters(Tcl_Interp *interp)
{
  int    mesh_max,                   mesh     = -1, tmp_mesh;
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL=0.0;
  int    cao_min, cao_max,           cao      = -1, tmp_cao;

  double                             alpha_L  = -1, tmp_alpha_L=0.0;
  double                             accuracy = -1, tmp_accuracy=0.0;
  double                            time_best=1e20, tmp_time;
  char
    b1[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b2[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b3[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 17];
 
  P3M_TRACE(fprintf(stderr,"%d: dipolar P3M_adaptive_tune_parameters\n",this_node));

  if (skin == -1) {
    Tcl_AppendResult(interp, "dipolar p3m cannot be tuned, since the skin is not yet set", (char *) NULL);
    return TCL_ERROR;
  }

  /* preparation */
  mpi_bcast_event(P3M_COUNT_DIPOLES);

  /* Print Status */
  sprintf(b1,"%.5e",p3m.Daccuracy);
  Tcl_AppendResult(interp, "dipolar P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);

  sprintf(b2,"%d",p3m_sum_dip_part); 

  Tcl_PrintDouble(interp, p3m_sum_mu2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);

  if (p3m_sum_dip_part == 0) {
    Tcl_AppendResult(interp, "no dipolar particles in the system, cannot tune dipolar P3M", (char *) NULL);
    return (TCL_ERROR);
  }
  
  
  /* parameter ranges */
  if (p3m.Dmesh[0] == 0 ) {
    double expo;
    expo = log(pow((double)p3m_sum_dip_part,(1.0/3.0)))/log(2.0);  

    tmp_mesh = (int)(pow(2.0,(double)((int)expo))+0.1);
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    mesh_max = tmp_mesh * 256;
    /* avoid using more than 1 GB of FFT arrays (per default, see config.h) */
    if (mesh_max > P3M_MAX_MESH)
      mesh_max = P3M_MAX_MESH;
  }
  else {
    sprintf(b1, "%d", p3m.Dmesh[0]);
    Tcl_AppendResult(interp, "fixed mesh ", b1, "\n", (char *)NULL);
    tmp_mesh = mesh_max = p3m.Dmesh[0];
  }

  if(p3m.Dr_cut_iL == 0.0) {
    r_cut_iL_min = 0;
    r_cut_iL_max = min_local_box_l/2 - skin;
    r_cut_iL_min *= box_l_i[0];
    r_cut_iL_max *= box_l_i[0];
  }
  else {
    sprintf(b1, "%f", p3m.Dr_cut_iL);
    Tcl_AppendResult(interp, "fixed r_cut_iL ", b1, "\n", (char *)NULL);
    r_cut_iL_min = r_cut_iL_max = p3m.Dr_cut_iL;
  }

  if(p3m.Dcao == 0) {
    cao_min = 1;
    cao_max = 7;
    cao = 3;
  }
  else {
    sprintf(b1, "%d", p3m.Dcao);
    Tcl_AppendResult(interp, "fixed cao ", b1, "\n", (char *)NULL);
    cao_min = cao_max = cao = p3m.Dcao;
  }

  Tcl_AppendResult(interp, "Dmesh cao Dr_cut_iL   Dalpha_L     Derr         Drs_err    Dks_err    time [ms]\n", (char *) NULL);

  /* mesh loop */
  for (;tmp_mesh <= mesh_max; tmp_mesh *= 2) {
    tmp_cao = cao;
    tmp_time = Dp3m_m_time(interp, tmp_mesh,
			  cao_min, cao_max, &tmp_cao,
			  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL,
			  &tmp_alpha_L, &tmp_accuracy);
    /* some error occured during the tuning force evaluation */
    if (tmp_time == -1) return TCL_ERROR;
    /* this mesh does not work at all */
    if (tmp_time < 0) continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      time_best = tmp_time;
      mesh      = tmp_mesh;
      cao       = tmp_cao;
      r_cut_iL  = tmp_r_cut_iL;
      alpha_L   = tmp_alpha_L;
      accuracy  = tmp_accuracy;
    }
    /* no hope of further optimisation */
    else if (tmp_time > time_best + P3M_TIME_GRAN)
      break;
  }
  
  P3M_TRACE(fprintf(stderr,"finshed tuning\n"));
  if(time_best == 1e20) {
    Tcl_AppendResult(interp, "failed to tune dipolar P3M parameters to required accuracy", (char *) NULL);
    return (TCL_ERROR);
  }

  /* set tuned p3m parameters */
  p3m.Dr_cut_iL = r_cut_iL;
  p3m.Dmesh[0]  = p3m.Dmesh[1] = p3m.Dmesh[2] = mesh;
  p3m.Dcao      = cao;
  p3m.Dalpha_L  = alpha_L;
  p3m.Daccuracy = accuracy;
  P3M_scaleby_box_l_dipoles();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  Tcl_AppendResult(interp, "\nresulting parameters:\n", (char *) NULL);
  sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
  Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",alpha_L); sprintf(b3,"%.5e",accuracy);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
  sprintf(b3,"                 %-8d",(int)time_best);
  Tcl_AppendResult(interp, b3, (char *) NULL);
  return (TCL_OK);
}


/*****************************************************************************/

void P3M_count_magnetic_particles()
{  
  Cell *cell;
  Particle *part;
  int i,c,np;
  double node_sums[2], tot_sums[2];

  for(i=0;i<2;i++)
    { node_sums[i]=0.0; tot_sums[i]=0.0;}

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
     if( part[i].p.dipm != 0.0 ) {
	node_sums[0] += SQR(part[i].r.dip[0])
	               +SQR(part[i].r.dip[1])
		       +SQR(part[i].r.dip[2]);
	node_sums[1] += 1.0;	       
      }
    }
  }
  
  MPI_Allreduce(node_sums, tot_sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  p3m_sum_mu2 = tot_sums[0];
  p3m_sum_dip_part    = (int)(tot_sums[1]+0.1);  
}



/*****************************************************************************/

/* Following are the two functions for computing the error in dipolar P3M and 
   tune the parameters to minimize the time with the desired accuracy.
   
   
   This functions are called by the functions: Dget_accuracy() and
   DP3M_tune_parameters.
  
*/


/*****************************************************************************/



   
double P3M_DIPOLAR_k_space_error(double box_size, double prefac, int mesh, int cao, int n_c_part, double sum_q2, double alpha_L)
{
  int  nx, ny, nz;
  double he_q = 0.0, mesh_i = 1./mesh, alpha_L_i = 1./alpha_L;
  double alias1, alias2, n2, cs;


 
  for (nx=-mesh/2; nx<mesh/2; nx++)
    for (ny=-mesh/2; ny<mesh/2; ny++)
      for (nz=-mesh/2; nz<mesh/2; nz++)
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  n2 = SQR(nx) + SQR(ny) + SQR(nz);
	  cs = analytic_cotangent_sum(nx,mesh_i,cao)*
 	       analytic_cotangent_sum(ny,mesh_i,cao)*
	       analytic_cotangent_sum(nz,mesh_i,cao);
	  P3M_DIPOLAR_tune_aliasing_sums(nx,ny,nz,mesh,mesh_i,cao,alpha_L_i,&alias1,&alias2);
	  he_q += (alias1  -  SQR(alias2/cs) / (n2*n2*n2));
	}

  return 8.*PI*PI/3.*sum_q2*sqrt(he_q/(double)n_c_part) / (box_size*box_size*box_size*box_size);
}
   

void P3M_DIPOLAR_tune_aliasing_sums(int nx, int ny, int nz,  int mesh, double mesh_i, int cao, double alpha_L_i,  double *alias1, double *alias2)
{
  int    mx,my,mz;
  double nmx,nmy,nmz;
  double fnmx,fnmy,fnmz;

  double ex,ex2,nm2,U2,factor1;

  factor1 = SQR(PI*alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (mx=-P3M_BRILLOUIN; mx<=P3M_BRILLOUIN; mx++) {
    fnmx = mesh_i * (nmx = nx + mx*mesh);
    for (my=-P3M_BRILLOUIN; my<=P3M_BRILLOUIN; my++) {
      fnmy = mesh_i * (nmy = ny + my*mesh);
      for (mz=-P3M_BRILLOUIN; mz<=P3M_BRILLOUIN; mz++) {
	fnmz = mesh_i * (nmz = nz + mz*mesh);
	
	nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	ex2 = SQR( ex = exp(-factor1*nm2) );
	
	U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex2 * nm2;
	*alias2 += U2 * ex * pow((nx*nmx + ny*nmy + nz*nmz),3.) / nm2;
     }
    }
  }
}





//----------------------------------------------------------
//  Function used to calculate the value of the errors
//  for the REAL part of the force in terms of the Spliting parameter alpha of Ewald
//  based on the formulas 33, the paper of Zuowei-HolmJCP, 115,6351,(2001).

//  Please, notice than in this more refined approach we don't use
//  formulas 37, but 33 which mantains all the powers in alpha

//------------------------------------------------------------



   
 double P3M_DIPOLAR_real_space_error(double box_size, double prefac, double r_cut_iL,  int n_c_part, double sum_q2, double alpha_L)
{
  double d_error_f,d_cc,d_dc,d_bc,d_rcut2,d_con;
  double d_a2,d_c,d_RCUT;
 
   
   
   
  d_RCUT=r_cut_iL*box_size;
  d_rcut2=d_RCUT*d_RCUT;
  
  
 d_a2=alpha_L*alpha_L/(box_size*box_size);

 d_c =sum_q2*exp(-d_a2*d_RCUT*d_RCUT);

 d_cc=4.0*d_a2*d_a2*d_rcut2*d_rcut2+6.0*d_a2*d_rcut2+3.0;


 d_dc=8.0*d_a2*d_a2*d_a2*d_rcut2*d_rcut2*d_rcut2+20.0*d_a2*d_a2*d_rcut2*d_rcut2 \
      +30*d_a2*d_rcut2+15.0;

 d_bc=2.0*d_a2*d_rcut2 +1.0;


 d_con=1.0/sqrt(box_size*box_size*box_size*d_a2*d_a2*d_rcut2*d_rcut2*d_rcut2*d_rcut2*d_RCUT*(double)n_c_part);


 d_error_f=d_c*d_con*sqrt((13./6.)*d_cc*d_cc+(2./15.)*d_dc*d_dc-(13./15.)*d_cc*d_dc);
 
 
  return d_error_f;
}

  

/*****************************************************************************/
  
    
// Using bisection find the root of a function "func-tuned_accuracy/sqrt(2.)" known to lie
//between x1 and x2. The root, returned as rtbis, will be refined
//until its accuracy is +-xacc.

double JJ_rtbisection( double box_size, double prefac, double r_cut_iL,  int n_c_part, double sum_q2,  double x1, double x2, double xacc, double tuned_accuracy)
{
 int j;
 double dx,f,fmid,xmid,rtb,constant,JJ_RTBIS_MAX=40;
 
 constant=tuned_accuracy/sqrt(2.);
 
 
 f=P3M_DIPOLAR_real_space_error(box_size,prefac,r_cut_iL,n_c_part,sum_q2,       x1)-constant;
 fmid=P3M_DIPOLAR_real_space_error(box_size,prefac,r_cut_iL,n_c_part,sum_q2,    x2)-constant;
 if(f*fmid >=0.0) fprintf(stderr,"Root must be bracketed for bisection in JJ_rtbisection \n");
 rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);  // Orient the search dx, and set rtb to x1 or x2 ...
 for (j=1;j<=JJ_RTBIS_MAX;j++){
   fmid=P3M_DIPOLAR_real_space_error(box_size,prefac,r_cut_iL,n_c_part,sum_q2,  xmid=rtb+(dx *= 0.5))-constant;
   if(fmid<=0.0) rtb=xmid;
   if(fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  fprintf(stderr,"Too many bisections in JJ_rtbissection \n");
  return -9999999.9999;  

}       
 
 

/************************************************************/


void Dcalc_lm_ld_pos() {
  int i; 
   for(i=0;i<3;i++) {
    Dlm.ld_pos[i] = (Dlm.ld_ind[i]+ p3m.Dmesh_off[i])*p3m.Da[i];
  }
}



/************************************************************/


void DP3M_init_a_ai_cao_cut() {
  int i; 
   for(i=0;i<3;i++) {
    p3m.Dai[i]      = (double)p3m.Dmesh[i]/box_l[i]; 
    p3m.Da[i]       = 1.0/p3m.Dai[i];
    p3m.Dcao_cut[i] = 0.5*p3m.Da[i]*p3m.Dcao;
  }
}



/************************************************************/

void Dcalc_local_ca_mesh() {
  int i;
  int ind[3];
  /* total skin size */
  double full_skin[3];
  
   for(i=0;i<3;i++)
    full_skin[i]= p3m.Dcao_cut[i]+skin+p3m.Dadditional_mesh[i];

  /* inner left down grid point (global index) */
  for(i=0;i<3;i++) Dlm.in_ld[i] = (int)ceil(my_left[i]*p3m.Dai[i]-p3m.Dmesh_off[i]);
  /* inner up right grid point (global index) */
  for(i=0;i<3;i++) Dlm.in_ur[i] = (int)floor(my_right[i]*p3m.Dai[i]-p3m.Dmesh_off[i]);
  
  /* correct roundof errors at boundary */
  for(i=0;i<3;i++) {
    if((my_right[i]*p3m.Dai[i]-p3m.Dmesh_off[i])-Dlm.in_ur[i]<ROUND_ERROR_PREC) Dlm.in_ur[i]--;
    if(1.0+(my_left[i]*p3m.Dai[i]-p3m.Dmesh_off[i])-Dlm.in_ld[i]<ROUND_ERROR_PREC) Dlm.in_ld[i]--;
  }
  /* inner grid dimensions */
  for(i=0;i<3;i++) Dlm.inner[i] = Dlm.in_ur[i] - Dlm.in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for(i=0;i<3;i++) 
    Dlm.ld_ind[i]=(int)ceil((my_left[i]-full_skin[i])*p3m.Dai[i]-p3m.Dmesh_off[i]);
  /* spacial position of left down mesh point */
  Dcalc_lm_ld_pos();
  /* left down margin */
  for(i=0;i<3;i++) Dlm.margin[i*2] = Dlm.in_ld[i]-Dlm.ld_ind[i];
  /* up right grid point */
  for(i=0;i<3;i++) ind[i]=(int)floor((my_right[i]+full_skin[i])*p3m.Dai[i]-p3m.Dmesh_off[i]);
  /* correct roundof errors at up right boundary */
  for(i=0;i<3;i++)
    if(((my_right[i]+full_skin[i])*p3m.Dai[i]-p3m.Dmesh_off[i])-ind[i]==0) ind[i]--;
  /* up right margin */
  for(i=0;i<3;i++) Dlm.margin[(i*2)+1] = ind[i] - Dlm.in_ur[i];

  /* grid dimension */
  Dlm.size=1; 
  for(i=0;i<3;i++) {Dlm.dim[i] = ind[i] - Dlm.ld_ind[i] + 1; Dlm.size*=Dlm.dim[i];}
  /* reduce inner grid indices from global to local */
  for(i=0;i<3;i++) Dlm.in_ld[i] = Dlm.margin[i*2];
  for(i=0;i<3;i++) Dlm.in_ur[i] = Dlm.margin[i*2]+Dlm.inner[i];

  Dlm.q_2_off  = Dlm.dim[2] - p3m.Dcao;
  Dlm.q_21_off = Dlm.dim[2] * (Dlm.dim[1] - p3m.Dcao);
   
}

/*****************************************************************************/


int DP3M_sanity_checks_boxl() {
  char *errtxt;
  int i, ret = 0;
     for(i=0;i<3;i++) {
    /* check k-space cutoff */
    if(p3m.Dcao_cut[i] >= 0.5*box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      ERROR_SPRINTF(errtxt,"{039 dipolar P3M_init: k-space cutoff %f is larger than half of box dimension %f} ",p3m.Dcao_cut[i],box_l[i]);
      ret = 1;
    }
    if(p3m.Dcao_cut[i] >= local_box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      ERROR_SPRINTF(errtxt,"{040 dipolar P3M_init: k-space cutoff %f is larger than local box dimension %f} ",p3m.Dcao_cut[i],local_box_l[i]);
      ret = 1;
    }
  }
   return ret;
}

/*****************************************************************************/


int DP3M_sanity_checks()
{
  char *errtxt;
  int ret = 0;

  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{041 dipolar P3M requires periodicity 1 1 1} ");
    ret = 1;
  }
  
  if (n_nodes != 1) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{110 dipolar P3M does not run in parallel} ");
    ret = 1;
  }
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{042 dipolar P3M at present requires the domain decomposition cell system} ");
    ret = 1;
  }
  
  if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{043 dipolar P3M requires a cubic box} ");
    ret = 1;
  }

    if( (p3m.Dmesh[0] != p3m.Dmesh[1]) || (p3m.Dmesh[1] != p3m.Dmesh[2]) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{044 dipolar P3M requires a cubic mesh} ");
    ret = 1;
  }
  
  

  if (DP3M_sanity_checks_boxl()) ret = 1;

  if (skin == -1) {
    errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
    ERROR_SPRINTF(errtxt,"{047 dipolar P3M_init: skin is not yet set} ");
    ret = 1;
  }
  
  if( p3m.Dmesh[0] == 0) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{045 dipolar P3M_init: mesh size is not yet set} ");
    ret = 1;
  }
  if( p3m.Dcao == 0) {
    errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
    ERROR_SPRINTF(errtxt,"{046 dipolar P3M_init: cao is not yet set} ");
    ret = 1;
  }
  
  
  return ret;
}


/*****************************************************************************/


void Dcalc_send_mesh()
{
  int i,j,evenodd;
  int done[3]={0,0,0};
  MPI_Status status;
  /* send grids */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      /* left */
      Dsm.s_ld[i*2][j] = 0 + done[j]*Dlm.margin[j*2];
      if(j==i) Dsm.s_ur[i*2][j] = Dlm.margin[j*2]; 
      else     Dsm.s_ur[i*2][j] = Dlm.dim[j]-done[j]*Dlm.margin[(j*2)+1];
      /* right */
      if(j==i) Dsm.s_ld[(i*2)+1][j] = Dlm.in_ur[j];
      else     Dsm.s_ld[(i*2)+1][j] = 0 + done[j]*Dlm.margin[j*2];
      Dsm.s_ur[(i*2)+1][j] = Dlm.dim[j] - done[j]*Dlm.margin[(j*2)+1];
    }   
    done[i]=1;
  }
  Dsm.max=0;
  for(i=0;i<6;i++) {
    Dsm.s_size[i] = 1;
    for(j=0;j<3;j++) {
      Dsm.s_dim[i][j] = Dsm.s_ur[i][j]-Dsm.s_ld[i][j];
      Dsm.s_size[i] *= Dsm.s_dim[i][j];
    }
    if(Dsm.s_size[i]>Dsm.max) Dsm.max=Dsm.s_size[i];
  }
  /* communication */
  for(i=0;i<6;i++) {
    if(i%2==0) j = i+1;
    else       j = i-1;
    if(node_neighbors[i] != this_node) {
      /* two step communication: first all even positions than all odd */
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[i/2]+evenodd)%2==0)
	  MPI_Send(&(Dlm.margin[i]), 1, MPI_INT, 
		   node_neighbors[i],REQ_P3M_INIT_D,MPI_COMM_WORLD);
	else
	  MPI_Recv(&(Dlm.r_margin[j]), 1, MPI_INT,
		   node_neighbors[j],REQ_P3M_INIT_D,MPI_COMM_WORLD,&status);    
      }
    }
    else {
      Dlm.r_margin[j] = Dlm.margin[i];
    }
  }
  /* recv grids */
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) {
      if(j==i) {
	Dsm.r_ld[ i*2   ][j] = Dsm.s_ld[ i*2   ][j] + Dlm.margin[2*j];
	Dsm.r_ur[ i*2   ][j] = Dsm.s_ur[ i*2   ][j] + Dlm.r_margin[2*j];
	Dsm.r_ld[(i*2)+1][j] = Dsm.s_ld[(i*2)+1][j] - Dlm.r_margin[(2*j)+1];
	Dsm.r_ur[(i*2)+1][j] = Dsm.s_ur[(i*2)+1][j] - Dlm.margin[(2*j)+1];
      }
      else {
	Dsm.r_ld[ i*2   ][j] = Dsm.s_ld[ i*2   ][j];
	Dsm.r_ur[ i*2   ][j] = Dsm.s_ur[ i*2   ][j];
	Dsm.r_ld[(i*2)+1][j] = Dsm.s_ld[(i*2)+1][j];
	Dsm.r_ur[(i*2)+1][j] = Dsm.s_ur[(i*2)+1][j];
      }
    }
  for(i=0;i<6;i++) {
    Dsm.r_size[i] = 1;
    for(j=0;j<3;j++) {
      Dsm.r_dim[i][j] = Dsm.r_ur[i][j]-Dsm.r_ld[i][j];
      Dsm.r_size[i] *= Dsm.r_dim[i][j];
    }
    if(Dsm.r_size[i]>Dsm.max) Dsm.max=Dsm.r_size[i];
  }
  
  
}



/************************************************/

void P3M_scaleby_box_l_dipoles() {
 
  p3m.Dr_cut = p3m.Dr_cut_iL* box_l[0];
  p3m.Dalpha = p3m.Dalpha_L * box_l_i[0];  
  DP3M_init_a_ai_cao_cut();
  Dcalc_lm_ld_pos();
  DP3M_sanity_checks_boxl();
  Dflag_constants_energy_dipolar=0; /* to ensure constants will be computed in case you want to calculate the energy */
}

/*****************************************************************************/
 
 
 /* fucntion to give the dipolar-P3M energy  correction -------*/
 void compute_constants_energy_dipolar() {
      double  Eself, Ukp3m;
      double volume;
 
       if(Dflag_constants_energy_dipolar==1) { 
 
 /*             double sumi1_value,sumi2_value, Eself, Ukp3m, Uk, Ur_cut;
                sumi1_value=JJ_sumi1(p3m.Dalpha_L);
            sumi2_value=JJ_sumi2(p3m.Dalpha_L);
             Eself=-1.*(*2*pow(p3m.Dalpha_L*box_l_i[0],3) * wupii/3.0);
	     Ukp3m=P3M_Average_dipolar_SelfEnergy(box_l[0],p3m.Dmesh[0]);
	     Uk=4.*PI*p3m_sum_mu2/6./pow(box_l[0],3)*sumi1_value;
	     Ur_cut=-p3m_sum_mu2*4.*pow(p3m.Dalpha_L*box_l_i[0],3)*wupii/6.*sumi2_value;	
	     Dipolar_energy_correction= Eself + Uk - Ukp3m + Ur_cut;
	     */
	     
	     Ukp3m=P3M_Average_dipolar_SelfEnergy(box_l[0],p3m.Dmesh[0]);
             Eself=-(2*pow(p3m.Dalpha_L*box_l_i[0],3) * wupii/3.0);
             volume=box_l[0]*box_l[1]*box_l[2];

       	     Dipolar_energy_correction= - p3m_sum_mu2*(Ukp3m+Eself+2.*PI/(3.*volume));
	     
	     /*fprintf(stderr,"p3m_sum_mu2=%lf\n",p3m_sum_mu2);
	     fprintf(stderr,"Ukp3m=%lf\n",Ukp3m);
	     fprintf(stderr,"Eself=%lf\n",Eself);
	     fprintf(stderr,"2.*PI/(3.*volume)=%lf\n",2.*PI/(3.*volume));*/
	 
	    Dflag_constants_energy_dipolar=2;
        }
  } 

/*****************************************************************************/


void   P3M_init_dipoles() {

  int n;

  if (coulomb.Dbjerrum == 0.0) {       
       if(coulomb.Dbjerrum == 0.0) {
           p3m.Dr_cut    = 0.0;
           p3m.Dr_cut_iL = 0.0;
          if(this_node==0) 
             P3M_TRACE(fprintf(stderr,"0: P3M_init_dipoles: dipolar Bjerrum length is zero.\n");
	   fprintf(stderr,"   Magnetostatics of dipoles switched off!\n"));
      }
  } else {  
    P3M_TRACE(fprintf(stderr,"%d: P3M_init_dipoles: \n",this_node));

    if (DP3M_sanity_checks()) return;

    P3M_TRACE(fprintf(stderr,"%d: P3M_init_dipoles: starting\n",this_node));

        P3M_TRACE(fprintf(stderr,"%d: mesh=%d, cao=%d, mesh_off=(%f,%f,%f)\n",this_node,p3m.Dmesh[0],p3m.Dcao,p3m.Dmesh_off[0],p3m.Dmesh_off[1],p3m.Dmesh_off[2]));
        p3m.Dcao3 = p3m.Dcao*p3m.Dcao*p3m.Dcao;


    /* initializes the (inverse) mesh constant p3m.Da (p3m.Dai) and the cutoff for charge assignment p3m.Dcao_cut */
    DP3M_init_a_ai_cao_cut();

    /* initialize ca fields to size CA_INCREMENT: Dca_frac and Dca_fmp */
    Dca_num = 0;
    if(Dca_num < CA_INCREMENT) {
      Dca_num = 0;
      Drealloc_ca_fields(CA_INCREMENT);
    }
 
    Dcalc_local_ca_mesh();

       Dcalc_send_mesh();
       P3M_TRACE(p3m_print_local_mesh(Dlm));
    
    /* DEBUG */
    for(n=0;n<n_nodes;n++) {
      /* MPI_Barrier(MPI_COMM_WORLD); */
         if(n==this_node) P3M_TRACE(p3m_print_send_mesh(Dsm));
    }
    
    if(Dsm.max != Dsend_recv_grid_size) {
      Dsend_recv_grid_size= Dsm.max;
      Dsend_grid = (double *) realloc(Dsend_grid, sizeof(double)*Dsm.max);
      Drecv_grid = (double *) realloc(Drecv_grid, sizeof(double)*Dsm.max);
    }
    
     if (p3m.Dinter > 0) interpolate_dipole_assignment_function();

     Dpos_shift = (double)((p3m.Dcao-1)/2) - (p3m.Dcao%2)/2.0;
     P3M_TRACE(fprintf(stderr,"%d: dipolar pos_shift = %f\n",this_node,Dpos_shift)); 
 
    /* FFT */
      P3M_TRACE(fprintf(stderr,"%d: Drs_mesh ADR=%p\n",this_node,Drs_mesh));
 
    Dca_mesh_size = Dfft_init(&Drs_mesh,Dlm.dim,Dlm.margin,&Dks_pnum);
    Dks_mesh = (double *) realloc(Dks_mesh, Dca_mesh_size*sizeof(double));

    for (n=0;n<3;n++)   
       Drs_mesh_dip[n] = (double *) realloc(Drs_mesh_dip[n], Dca_mesh_size*sizeof(double));

     P3M_TRACE(fprintf(stderr,"%d: Drs_mesh_dip[0] ADR=%p\n",this_node,Drs_mesh_dip[0]));
     P3M_TRACE(fprintf(stderr,"%d: Drs_mesh_dip[1] ADR=%p\n",this_node,Drs_mesh_dip[1]));
     P3M_TRACE(fprintf(stderr,"%d: Drs_mesh_dip[2] ADR=%p\n",this_node,Drs_mesh_dip[2]));
 
 
    /* k-space part: */
    
    Dcalc_differential_operator();

    Dcalc_influence_function_force();
    Dcalc_influence_function_energy();

    P3M_count_magnetic_particles();
    
    Dflag_constants_energy_dipolar=0; /* to ensure constants will be computed in case you want to calculate the energy */
   
    P3M_TRACE(fprintf(stderr,"%d: p3m initialized\n",this_node));
  }


}

/*****************************************************************************/


#endif  /* of defined(MAGNETOSTATICS) */

