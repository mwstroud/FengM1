#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _bg2pyr_reg(void);
extern void _cadyn_reg(void);
extern void _cal2CA3_reg(void);
extern void _cal2_reg(void);
extern void _ca_reg(void);
extern void _capoolCA3_reg(void);
extern void _capool_reg(void);
extern void _cas_reg(void);
extern void _cat_reg(void);
extern void _currentclamp_reg(void);
extern void _function_TMonitor_reg(void);
extern void _gap_reg(void);
extern void _Gfluct_new_exc_reg(void);
extern void _Gfluct_new_inh_reg(void);
extern void _hCA3_reg(void);
extern void _h_reg(void);
extern void _imCA3_reg(void);
extern void _im_reg(void);
extern void _interD2interD_SOMPV_STFD_new_reg(void);
extern void _interD2interD_STFD_new_reg(void);
extern void _interD2pyrD_SOM2P_STFD_new_reg(void);
extern void _interD2pyrD_STFD_new_reg(void);
extern void _kadist_reg(void);
extern void _kaprox_reg(void);
extern void _kca_reg(void);
extern void _kdrca1DA_reg(void);
extern void _kdrca1_reg(void);
extern void _kdrCA3_reg(void);
extern void _kdrinter_reg(void);
extern void _leakCA3_reg(void);
extern void _leakDA_reg(void);
extern void _leakinter_reg(void);
extern void _leak_reg(void);
extern void _na3DA_reg(void);
extern void _na3_reg(void);
extern void _nainter_reg(void);
extern void _napCA3_reg(void);
extern void _nap_reg(void);
extern void _natCA3_reg(void);
extern void _nat_reg(void);
extern void _pyrD2interD_P2SOM_STFD_reg(void);
extern void _pyrD2interD_STFD_reg(void);
extern void _pyrD2pyrD_STFD_new_reg(void);
extern void _sahpCA3_reg(void);
extern void _sahp_reg(void);
extern void _sahpNE_reg(void);
extern void _vecevent_reg(void);
extern void _xtra_imemrec_reg(void);
extern void _xtra_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," modfiles/bg2pyr.mod");
    fprintf(stderr," modfiles/cadyn.mod");
    fprintf(stderr," modfiles/cal2CA3.mod");
    fprintf(stderr," modfiles/cal2.mod");
    fprintf(stderr," modfiles/ca.mod");
    fprintf(stderr," modfiles/capoolCA3.mod");
    fprintf(stderr," modfiles/capool.mod");
    fprintf(stderr," modfiles/cas.mod");
    fprintf(stderr," modfiles/cat.mod");
    fprintf(stderr," modfiles/currentclamp.mod");
    fprintf(stderr," modfiles/function_TMonitor.mod");
    fprintf(stderr," modfiles/gap.mod");
    fprintf(stderr," modfiles/Gfluct_new_exc.mod");
    fprintf(stderr," modfiles/Gfluct_new_inh.mod");
    fprintf(stderr," modfiles/hCA3.mod");
    fprintf(stderr," modfiles/h.mod");
    fprintf(stderr," modfiles/imCA3.mod");
    fprintf(stderr," modfiles/im.mod");
    fprintf(stderr," modfiles/interD2interD_SOMPV_STFD_new.mod");
    fprintf(stderr," modfiles/interD2interD_STFD_new.mod");
    fprintf(stderr," modfiles/interD2pyrD_SOM2P_STFD_new.mod");
    fprintf(stderr," modfiles/interD2pyrD_STFD_new.mod");
    fprintf(stderr," modfiles/kadist.mod");
    fprintf(stderr," modfiles/kaprox.mod");
    fprintf(stderr," modfiles/kca.mod");
    fprintf(stderr," modfiles/kdrca1DA.mod");
    fprintf(stderr," modfiles/kdrca1.mod");
    fprintf(stderr," modfiles/kdrCA3.mod");
    fprintf(stderr," modfiles/kdrinter.mod");
    fprintf(stderr," modfiles/leakCA3.mod");
    fprintf(stderr," modfiles/leakDA.mod");
    fprintf(stderr," modfiles/leakinter.mod");
    fprintf(stderr," modfiles/leak.mod");
    fprintf(stderr," modfiles/na3DA.mod");
    fprintf(stderr," modfiles/na3.mod");
    fprintf(stderr," modfiles/nainter.mod");
    fprintf(stderr," modfiles/napCA3.mod");
    fprintf(stderr," modfiles/nap.mod");
    fprintf(stderr," modfiles/natCA3.mod");
    fprintf(stderr," modfiles/nat.mod");
    fprintf(stderr," modfiles/pyrD2interD_P2SOM_STFD.mod");
    fprintf(stderr," modfiles/pyrD2interD_STFD.mod");
    fprintf(stderr," modfiles/pyrD2pyrD_STFD_new.mod");
    fprintf(stderr," modfiles/sahpCA3.mod");
    fprintf(stderr," modfiles/sahp.mod");
    fprintf(stderr," modfiles/sahpNE.mod");
    fprintf(stderr," modfiles/vecevent.mod");
    fprintf(stderr," modfiles/xtra_imemrec.mod");
    fprintf(stderr," modfiles/xtra.mod");
    fprintf(stderr, "\n");
  }
  _bg2pyr_reg();
  _cadyn_reg();
  _cal2CA3_reg();
  _cal2_reg();
  _ca_reg();
  _capoolCA3_reg();
  _capool_reg();
  _cas_reg();
  _cat_reg();
  _currentclamp_reg();
  _function_TMonitor_reg();
  _gap_reg();
  _Gfluct_new_exc_reg();
  _Gfluct_new_inh_reg();
  _hCA3_reg();
  _h_reg();
  _imCA3_reg();
  _im_reg();
  _interD2interD_SOMPV_STFD_new_reg();
  _interD2interD_STFD_new_reg();
  _interD2pyrD_SOM2P_STFD_new_reg();
  _interD2pyrD_STFD_new_reg();
  _kadist_reg();
  _kaprox_reg();
  _kca_reg();
  _kdrca1DA_reg();
  _kdrca1_reg();
  _kdrCA3_reg();
  _kdrinter_reg();
  _leakCA3_reg();
  _leakDA_reg();
  _leakinter_reg();
  _leak_reg();
  _na3DA_reg();
  _na3_reg();
  _nainter_reg();
  _napCA3_reg();
  _nap_reg();
  _natCA3_reg();
  _nat_reg();
  _pyrD2interD_P2SOM_STFD_reg();
  _pyrD2interD_STFD_reg();
  _pyrD2pyrD_STFD_new_reg();
  _sahpCA3_reg();
  _sahp_reg();
  _sahpNE_reg();
  _vecevent_reg();
  _xtra_imemrec_reg();
  _xtra_reg();
}
