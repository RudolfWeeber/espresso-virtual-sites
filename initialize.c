#include "initialize.h"
#include "interaction_data.h"
#include "binary_file.h"
#include "integrate.h"

int initialize(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
  */


  /*
    installation of tcl commands
  */

  /* in integrate.c */
  Tcl_CreateCommand(interp, "integrate", integrate, 0, NULL);
  /* in global.c */
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "inter", inter, 0, NULL);
  /* in particle_data.c */
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
  /* in file binaryfile.c */
  Tcl_CreateCommand(interp, "writemd", writemd, 0, NULL);
  Tcl_CreateCommand(interp, "readmd", readmd, 0, NULL);

  return TCL_OK;
}
