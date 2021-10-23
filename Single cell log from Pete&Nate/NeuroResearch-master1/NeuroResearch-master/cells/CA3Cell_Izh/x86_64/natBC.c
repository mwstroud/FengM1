/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__natBC
#define _nrn_initial _nrn_initial__natBC
#define nrn_cur _nrn_cur__natBC
#define _nrn_current _nrn_current__natBC
#define nrn_jacob _nrn_jacob__natBC
#define nrn_state _nrn_state__natBC
#define _net_receive _net_receive__natBC 
#define states states__natBC 
#define trates trates__natBC 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gbar _p[0]
#define km _p[1]
#define kh _p[2]
#define vhalfm _p[3]
#define vhalfh _p[4]
#define i _p[5]
#define gna _p[6]
#define minf _p[7]
#define hinf _p[8]
#define mtau _p[9]
#define htau _p[10]
#define m _p[11]
#define h _p[12]
#define ena _p[13]
#define ina _p[14]
#define tha1 _p[15]
#define Dm _p[16]
#define Dh _p[17]
#define v _p[18]
#define _g _p[19]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_trap0(void);
 static void _hoc_trates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_natBC", _hoc_setdata,
 "trap0_natBC", _hoc_trap0,
 "trates_natBC", _hoc_trates,
 0, 0
};
#define trap0 trap0_natBC
 extern double trap0( _threadargsprotocomma_ double , double , double , double );
 #define _zmexp _thread[0]._pval[0]
 #define _zhexp _thread[0]._pval[1]
 /* declare global and static user variables */
#define Rd Rd_natBC
 double Rd = 0.029;
#define Rg Rg_natBC
 double Rg = 0.01;
#define Rb Rb_natBC
 double Rb = 0.5;
#define Ra Ra_natBC
 double Ra = 0.4;
#define ar2 ar2_natBC
 double ar2 = 1;
#define hmin hmin_natBC
 double hmin = 0.5;
#define mmin mmin_natBC
 double mmin = 0.02;
#define qinf qinf_natBC
 double qinf = 4;
#define qq qq_natBC
 double qq = 10;
#define q10 q10_natBC
 double q10 = 2;
#define qg qg_natBC
 double qg = 1.5;
#define qd qd_natBC
 double qd = 1.5;
#define qa qa_natBC
 double qa = 10;
#define thinf thinf_natBC
 double thinf = -50;
#define tq tq_natBC
 double tq = -55;
#define thi2 thi2_natBC
 double thi2 = -35;
#define thi1 thi1_natBC
 double thi1 = -40;
#define tha tha_natBC
 double tha = -30;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tha_natBC", "mV",
 "qa_natBC", "mV",
 "Ra_natBC", "/ms",
 "Rb_natBC", "/ms",
 "thi1_natBC", "mV",
 "thi2_natBC", "mV",
 "qd_natBC", "mV",
 "qg_natBC", "mV",
 "Rg_natBC", "/ms",
 "Rd_natBC", "/ms",
 "qq_natBC", "mV",
 "tq_natBC", "mV",
 "thinf_natBC", "mV",
 "qinf_natBC", "mV",
 "ar2_natBC", "1",
 "gbar_natBC", "mho/cm2",
 "km_natBC", "ms",
 "kh_natBC", "ms",
 "vhalfm_natBC", "mV",
 "vhalfh_natBC", "mV",
 "i_natBC", "mA/cm2",
 "gna_natBC", "mho/cm2",
 "mtau_natBC", "ms",
 "htau_natBC", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tha_natBC", &tha_natBC,
 "qa_natBC", &qa_natBC,
 "Ra_natBC", &Ra_natBC,
 "Rb_natBC", &Rb_natBC,
 "thi1_natBC", &thi1_natBC,
 "thi2_natBC", &thi2_natBC,
 "qd_natBC", &qd_natBC,
 "qg_natBC", &qg_natBC,
 "mmin_natBC", &mmin_natBC,
 "hmin_natBC", &hmin_natBC,
 "q10_natBC", &q10_natBC,
 "Rg_natBC", &Rg_natBC,
 "Rd_natBC", &Rd_natBC,
 "qq_natBC", &qq_natBC,
 "tq_natBC", &tq_natBC,
 "thinf_natBC", &thinf_natBC,
 "qinf_natBC", &qinf_natBC,
 "ar2_natBC", &ar2_natBC,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"natBC",
 "gbar_natBC",
 "km_natBC",
 "kh_natBC",
 "vhalfm_natBC",
 "vhalfh_natBC",
 0,
 "i_natBC",
 "gna_natBC",
 "minf_natBC",
 "hinf_natBC",
 "mtau_natBC",
 "htau_natBC",
 0,
 "m_natBC",
 "h_natBC",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 20, _prop);
 	/*initialize range parameters*/
 	gbar = 0.01;
 	km = 7.2;
 	kh = 4;
 	vhalfm = 38.43;
 	vhalfh = 50;
 	_prop->param = _p;
 	_prop->param_size = 20;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _natBC_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 20, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 natBC /home/pbcanfield/Desktop/NeuroResearch/CA3Cell_Izh/x86_64/natBC.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 /*Top LOCAL _zmexp , _zhexp */
static int _reset;
static char *modelname = "nat";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int trates(_threadargsprotocomma_ double, double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   trates ( _threadargscomma_ v , ar2 ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 trates ( _threadargscomma_ v , ar2 ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   trates ( _threadargscomma_ v , ar2 ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 
static int  trates ( _threadargsprotocomma_ double _lvm , double _la2 ) {
   double _la , _lb , _lqt ;
 _lqt = 1.6245 ;
   tha1 = tha ;
   _la = trap0 ( _threadargscomma_ _lvm , tha1 , Ra , qa ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - tha1 , Rb , qa ) ;
   mtau = ( 1.0 / ( _la + _lb ) / _lqt ) ;
   if ( mtau < mmin ) {
     mtau = mmin ;
     }
   if ( v < - 60.0 ) {
     minf = 0.0 ;
     }
   else {
     minf = 1.0 / ( 1.0 + exp ( ( - v - vhalfm ) / km ) ) ;
     }
   _la = trap0 ( _threadargscomma_ _lvm , thi1 , Rd , qd ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - thi2 , Rg , qg ) ;
   htau = 1.0 / ( _la + _lb ) / _lqt ;
   if ( htau < hmin ) {
     htau = hmin ;
     }
   hinf = 1.0 / ( 1.0 + exp ( ( v + vhalfh ) / kh ) ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 trates ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double trap0 ( _threadargsprotocomma_ double _lv , double _lth , double _la , double _lq ) {
   double _ltrap0;
 if ( fabs ( _lv - _lth ) > 1e-6 ) {
     _ltrap0 = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
     }
   else {
     _ltrap0 = _la * _lq ;
     }
   
return _ltrap0;
 }
 
static void _hoc_trap0(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  trap0 ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[0]._pval = (double*)ecalloc(2, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[0]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
   m = minf ;
   h = hinf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gbar * m * m * m * h ;
   ina = gna * ( v - ena ) ;
   i = ina ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/pbcanfield/Desktop/NeuroResearch/CA3Cell_Izh/modfiles/natBC.mod";
static const char* nmodl_file_text = 
  "TITLE nat\n"
  ": Na current \n"
  ": from Jeff M.\n"
  ":  ---------- modified -------M.Migliore may97\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX natBC\n"
  "	USEION na READ ena WRITE ina\n"
  "	:RANGE  , i :, ar2\n"
  "	RANGE gbar, gna, i, minf, hinf, mtau, htau, km, kh, vhalfh, vhalfm : , qinf, thinf\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	\n"
  "	gbar = 0.010   	(mho/cm2)	\n"
  "								\n"
  "	tha  =  -30	(mV)		: v 1/2 for act	\n"
  "	qa   = 10	(mV)		: act slope (4.5)		\n"
  "	Ra   = .4	(/ms)		: open (v)		\n"
  "	Rb   = 0.5 	(/ms)		: close (v) (0.124)	\n"
  "\n"
  "	thi1  = -40	(mV)		: v 1/2 for inact 	\n"
  "	thi2  = -35 	(mV)	: v 1/2 for inact 	\n"
  "	qd   = 1.5		(mV)    : inact tau slope\n"
  "	qg   = 1.5      (mV)\n"
  "	mmin=0.02	\n"
  "	hmin=0.5			\n"
  "	q10=2\n"
  "	Rg   = 0.01 	(/ms)		: inact recov (v) 	\n"
  "	Rd   = .029 	(/ms)		: inact (v)	\n"
  "	qq   = 10        (mV)\n"
  "	tq   = -55      (mV)\n"
  "\n"
  "	thinf  = -50 	(mV)		: inact inf slope	\n"
  "	qinf  = 4 	(mV)		: inact inf slope \n"
  "\n"
  "    ar2=1		(1)		: 1=no inact., 0=max inact.\n"
  "	ena		(mV)            : must be explicitly def. in hoc\n"
  "	celsius\n"
  "	v 		(mV)\n"
  "\n"
  "	km = 7.2 	(ms)\n"
  "	kh = 4 		(ms)\n"
  "\n"
  "	vhalfm = 38.43  (mV)\n"
  "	vhalfh = 50	(mV)\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "	ina 		(mA/cm2)\n"
  "	i    		(mA/cm2)\n"
  "	gna		(mho/cm2)\n"
  "	minf 		hinf 		\n"
  "	mtau (ms)	htau (ms) 	\n"
  "	tha1	\n"
  "}\n"
  " \n"
  "\n"
  "STATE { m h}\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gna = gbar*m*m*m*h\n"
  "	ina = gna * (v - ena)\n"
  "	i = ina\n"
  "} \n"
  "\n"
  "INITIAL {\n"
  "	trates(v,ar2)\n"
  "	m=minf  \n"
  "	h=hinf\n"
  "}\n"
  "\n"
  "\n"
  "LOCAL mexp, hexp\n"
  "\n"
  "DERIVATIVE states {   \n"
  "        trates(v,ar2)      \n"
  "        m' = (minf-m)/mtau\n"
  "        h' = (hinf-h)/htau\n"
  "}\n"
  "\n"
  "PROCEDURE trates(vm,a2) {  \n"
  "        LOCAL  a, b, qt\n"
  "		qt = 1.6245\n"
  "		tha1 = tha \n"
  "	a = trap0(vm,tha1,Ra,qa)\n"
  "	b = trap0(-vm,-tha1,Rb,qa)\n"
  "	mtau = (1/(a+b)/qt)\n"
  "        if (mtau<mmin) {mtau=mmin}\n"
  "	if (v < -60 ) {\n"
  "	minf = 0\n"
  "	} else{\n"
  "	minf  = 1 / ( 1 + exp( ( - v - vhalfm ) / km ) )\n"
  "	}\n"
  "	a = trap0(vm,thi1,Rd,qd)\n"
  "	b = trap0(-vm,-thi2,Rg,qg)\n"
  "	htau =  1/(a+b)/qt\n"
  "        \n"
  "	if (htau<hmin) {htau=hmin}\n"
  "	hinf  = 1 / ( 1 + exp( ( v + vhalfh ) / kh ) )\n"
  "}\n"
  "\n"
  "FUNCTION trap0(v,th,a,q) {\n"
  "	if (fabs(v-th) > 1e-6) {\n"
  "	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))\n"
  "	} else {\n"
  "	        trap0 = a * q\n"
  " 	}\n"
  "}	\n"
  ;
#endif
