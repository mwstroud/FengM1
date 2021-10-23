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
 
#define nrn_init _nrn_init__nca
#define _nrn_initial _nrn_initial__nca
#define nrn_cur _nrn_cur__nca
#define _nrn_current _nrn_current__nca
#define nrn_jacob _nrn_jacob__nca
#define nrn_state _nrn_state__nca
#define _net_receive _net_receive__nca 
#define rate rate__nca 
#define states states__nca 
 
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
#define gcabar _p[0]
#define uinf _p[1]
#define zinf _p[2]
#define utau _p[3]
#define ztau _p[4]
#define gca _p[5]
#define u _p[6]
#define z _p[7]
#define Du _p[8]
#define Dz _p[9]
#define eca _p[10]
#define ica _p[11]
#define v _p[12]
#define _g _p[13]
#define _ion_eca	*_ppvar[0]._pval
#define _ion_ica	*_ppvar[1]._pval
#define _ion_dicadv	*_ppvar[2]._pval
 
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
 /* declaration of user functions */
 static void _hoc_rate(void);
 static void _hoc_ubet(void);
 static void _hoc_ualf(void);
 static void _hoc_zbet(void);
 static void _hoc_zalf(void);
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
 "setdata_nca", _hoc_setdata,
 "rate_nca", _hoc_rate,
 "ubet_nca", _hoc_ubet,
 "ualf_nca", _hoc_ualf,
 "zbet_nca", _hoc_zbet,
 "zalf_nca", _hoc_zalf,
 0, 0
};
#define ubet ubet_nca
#define ualf ualf_nca
#define zbet zbet_nca
#define zalf zalf_nca
 extern double ubet( _threadargsprotocomma_ double );
 extern double ualf( _threadargsprotocomma_ double );
 extern double zbet( _threadargsprotocomma_ double );
 extern double zalf( _threadargsprotocomma_ double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gcabar_nca", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gcabar_nca", "siemens/cm2",
 "utau_nca", "ms",
 "ztau_nca", "ms",
 "gca_nca", "siemens/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double u0 = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
"nca",
 "gcabar_nca",
 0,
 "uinf_nca",
 "zinf_nca",
 "utau_nca",
 "ztau_nca",
 "gca_nca",
 0,
 "u_nca",
 "z_nca",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gcabar = 0.0001;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cell_nca_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 nca /home/pbcanfield/Desktop/NeuroResearch/CA3Cell_Izh/x86_64/cell_nca.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rate ( _threadargscomma_ v ) ;
   Du = ( uinf - u ) / utau ;
   Dz = ( zinf - z ) / ztau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rate ( _threadargscomma_ v ) ;
 Du = Du  / (1. - dt*( ( ( ( - 1.0 ) ) ) / utau )) ;
 Dz = Dz  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ztau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rate ( _threadargscomma_ v ) ;
    u = u + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / utau)))*(- ( ( ( uinf ) ) / utau ) / ( ( ( ( - 1.0 ) ) ) / utau ) - u) ;
    z = z + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ztau)))*(- ( ( ( zinf ) ) / ztau ) / ( ( ( ( - 1.0 ) ) ) / ztau ) - z) ;
   }
  return 0;
}
 
double ualf ( _threadargsprotocomma_ double _lv ) {
   double _lualf;
 double _lva ;
 _lva = 19.98 - _lv ;
   if ( fabs ( _lva ) < 1e-04 ) {
     _lva = _lva + 0.0001 ;
     }
   _lualf = 0.19 * _lva / ( exp ( _lva / 10.0 ) - 1.0 ) ;
   
return _lualf;
 }
 
static void _hoc_ualf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ualf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double ubet ( _threadargsprotocomma_ double _lv ) {
   double _lubet;
 double _lvb ;
 _lvb = _lv ;
   if ( fabs ( _lvb ) < 1e-04 ) {
     _lvb = _lvb + 0.0001 ;
     }
   _lubet = 0.046 * ( exp ( - _lvb / 20.73 ) ) ;
   
return _lubet;
 }
 
static void _hoc_ubet(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ubet ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double zalf ( _threadargsprotocomma_ double _lv ) {
   double _lzalf;
 double _lva ;
 _lva = _lv ;
   if ( fabs ( _lva ) < 1e-04 ) {
     _lva = _lva + 0.0001 ;
     }
   _lzalf = 1.6e-4 * exp ( - _lva / 48.4 ) ;
   
return _lzalf;
 }
 
static void _hoc_zalf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  zalf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double zbet ( _threadargsprotocomma_ double _lv ) {
   double _lzbet;
 double _lvb ;
 _lvb = 39.0 - _lv ;
   if ( fabs ( _lvb ) < 1e-04 ) {
     _lvb = _lvb + 0.0001 ;
     }
   _lzbet = 1.0 / ( exp ( _lvb / 10.0 ) + 1.0 ) ;
   
return _lzbet;
 }
 
static void _hoc_zbet(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  zbet ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rate ( _threadargsprotocomma_ double _lv ) {
   double _lusum , _lzsum , _lua , _lub , _lza , _lzb ;
 _lua = ualf ( _threadargscomma_ _lv ) ;
   _lub = ubet ( _threadargscomma_ _lv ) ;
   _lza = zalf ( _threadargscomma_ _lv ) ;
   _lzb = zbet ( _threadargscomma_ _lv ) ;
   _lusum = _lua + _lub ;
   uinf = _lua / _lusum ;
   utau = 1.0 / ( _lusum ) ;
   _lzsum = _lza + _lzb ;
   zinf = _lza / _lzsum ;
   ztau = 1.0 / ( _lzsum ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rate ( _p, _ppvar, _thread, _nt, *getarg(1) );
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
  eca = _ion_eca;
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
  eca = _ion_eca;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  u = u0;
  z = z0;
 {
   rate ( _threadargscomma_ v ) ;
   u = uinf ;
   z = zinf ;
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
  eca = _ion_eca;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gca = gcabar * u * u * z ;
   ica = gca * ( v - eca ) ;
   }
 _current += ica;

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
  eca = _ion_eca;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  eca = _ion_eca;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(u) - _p;  _dlist1[0] = &(Du) - _p;
 _slist1[1] = &(z) - _p;  _dlist1[1] = &(Dz) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/pbcanfield/Desktop/NeuroResearch/CA3Cell_Izh/modfiles/cell_nca.mod";
static const char* nmodl_file_text = 
  ":N-type voltage activated Ca current\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX nca\n"
  "	USEION ca READ eca WRITE ica\n"
  "	RANGE gcabar, gca\n"
  "	RANGE uinf, zinf, utau, ztau \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gcabar = 0.0001 (siemens/cm2) <0,1e9>\n"
  "}\n"
  "\n"
  "STATE { u z }\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	eca (mV)\n"
  "	ica (mA/cm2)\n"
  "	uinf\n"
  "	zinf \n"
  "	utau (ms)\n"
  "	ztau (ms)\n"
  "	gca (siemens/cm2)\n"
  "}\n"
  "\n"
  "BREAKPOINT { \n"
  "	SOLVE states METHOD cnexp\n"
  "	gca = gcabar*u*u*z\n"
  "	ica = gca*(v-eca)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rate(v)\n"
  "	u = uinf\n"
  "	z = zinf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rate(v)\n"
  "	u' = (uinf-u)/utau\n"
  "	z' = (zinf-z)/ztau\n"
  "}\n"
  "FUNCTION ualf(v(mV)) {\n"
  "	LOCAL va \n"
  "	va = 19.98-v						:v+25\n"
  "	if (fabs(va)<1e-04){\n"
  "		va = va+0.0001\n"
  "	}\n"
  "	ualf = 0.19*va/(exp(va/10)-1)\n"
  "}\n"
  "\n"
  "FUNCTION ubet(v(mV)) {\n"
  "	LOCAL vb \n"
  "	vb = v\n"
  "	if (fabs(vb)<1e-04){\n"
  "		vb = vb+0.0001\n"
  "	}\n"
  "	\n"
  "	ubet = 0.046*(exp(-vb/20.73))\n"
  "}	\n"
  "\n"
  "FUNCTION zalf(v(mV)) {\n"
  "	LOCAL va \n"
  "	va = v\n"
  "	if (fabs(va)<1e-04){\n"
  "		va = va+0.0001\n"
  "	}\n"
  "	\n"
  "	zalf = 1.6e-4*exp(-va/48.4)\n"
  "}\n"
  "\n"
  "FUNCTION zbet(v(mV)) {\n"
  "	LOCAL vb \n"
  "	vb = 39-v\n"
  "	if (fabs(vb)<1e-04){\n"
  "		vb = vb+0.0001\n"
  "	}\n"
  "	\n"
  "	zbet = 1/(exp(vb/10)+1)\n"
  "}\n"
  "PROCEDURE rate(v(mV)) { LOCAL usum, zsum, ua, ub, za, zb\n"
  "\n"
  "	ua = ualf(v) ub = ubet(v) za = zalf(v) zb = zbet(v)\n"
  "	\n"
  "	usum = ua+ub\n"
  "	uinf = ua/usum\n"
  "	utau = 1/(usum)\n"
  "	\n"
  "	zsum = za+zb\n"
  "	zinf = za/zsum\n"
  "	ztau = 1/(zsum)\n"
  "}\n"
  ;
#endif
