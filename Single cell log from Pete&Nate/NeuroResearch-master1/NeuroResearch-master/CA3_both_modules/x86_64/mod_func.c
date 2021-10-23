#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _cal2_reg(void);
extern void _capool_reg(void);
extern void _cas_reg(void);
extern void _cat_reg(void);
extern void _hCA3_reg(void);
extern void _im_reg(void);
extern void _kca_reg(void);
extern void _kdrCA3_reg(void);
extern void _leakCA3_reg(void);
extern void _nap_reg(void);
extern void _natCA3_reg(void);
extern void _sahp_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," cal2.mod");
    fprintf(stderr," capool.mod");
    fprintf(stderr," cas.mod");
    fprintf(stderr," cat.mod");
    fprintf(stderr," hCA3.mod");
    fprintf(stderr," im.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kdrCA3.mod");
    fprintf(stderr," leakCA3.mod");
    fprintf(stderr," nap.mod");
    fprintf(stderr," natCA3.mod");
    fprintf(stderr," sahp.mod");
    fprintf(stderr, "\n");
  }
  _cal2_reg();
  _capool_reg();
  _cas_reg();
  _cat_reg();
  _hCA3_reg();
  _im_reg();
  _kca_reg();
  _kdrCA3_reg();
  _leakCA3_reg();
  _nap_reg();
  _natCA3_reg();
  _sahp_reg();
}
