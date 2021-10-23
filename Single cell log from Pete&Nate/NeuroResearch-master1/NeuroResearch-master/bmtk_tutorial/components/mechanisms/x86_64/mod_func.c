#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _capool_reg(void);
extern void _cas_reg(void);
extern void _cat_reg(void);
extern void _hyper_reg(void);
extern void _inhsyn_reg(void);
extern void _ka_reg(void);
extern void _kca_reg(void);
extern void _kdr_reg(void);
extern void _leak_reg(void);
extern void _na_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," capool.mod");
    fprintf(stderr," cas.mod");
    fprintf(stderr," cat.mod");
    fprintf(stderr," hyper.mod");
    fprintf(stderr," inhsyn.mod");
    fprintf(stderr," ka.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kdr.mod");
    fprintf(stderr," leak.mod");
    fprintf(stderr," na.mod");
    fprintf(stderr, "\n");
  }
  _capool_reg();
  _cas_reg();
  _cat_reg();
  _hyper_reg();
  _inhsyn_reg();
  _ka_reg();
  _kca_reg();
  _kdr_reg();
  _leak_reg();
  _na_reg();
}
