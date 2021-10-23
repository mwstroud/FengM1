#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _kdrBC_reg(void);
extern void _leakBC_reg(void);
extern void _natBC_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," modfiles//kdrBC.mod");
    fprintf(stderr," modfiles//leakBC.mod");
    fprintf(stderr," modfiles//natBC.mod");
    fprintf(stderr, "\n");
  }
  _kdrBC_reg();
  _leakBC_reg();
  _natBC_reg();
}
