#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _OLM_IA_reg();
extern void _OLM_Ih_reg();
extern void _OLM_Lca_reg();
extern void _OLM_cat_reg();
extern void _OLM_ccanl_reg();
extern void _OLM_ichan2_reg();
extern void _OLM_sahp_reg();
extern void _bg2pyr_reg();
extern void _currentclamp_reg();
extern void _nap_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," OLM_IA.mod");
fprintf(stderr," OLM_Ih.mod");
fprintf(stderr," OLM_Lca.mod");
fprintf(stderr," OLM_cat.mod");
fprintf(stderr," OLM_ccanl.mod");
fprintf(stderr," OLM_ichan2.mod");
fprintf(stderr," OLM_sahp.mod");
fprintf(stderr," bg2pyr.mod");
fprintf(stderr," currentclamp.mod");
fprintf(stderr," nap.mod");
fprintf(stderr, "\n");
    }
_OLM_IA_reg();
_OLM_Ih_reg();
_OLM_Lca_reg();
_OLM_cat_reg();
_OLM_ccanl_reg();
_OLM_ichan2_reg();
_OLM_sahp_reg();
_bg2pyr_reg();
_currentclamp_reg();
_nap_reg();
}
