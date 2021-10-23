:N-type voltage activated Ca current

NEURON {
    SUFFIX nca
	USEION ca READ eca WRITE ica
	RANGE gcabar, gca
	RANGE uinf, zinf, utau, ztau 
}

UNITS {
    (mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcabar = 0.0001 (siemens/cm2) <0,1e9>
}

STATE { u z }

ASSIGNED {
	v (mV)
	eca (mV)
	ica (mA/cm2)
	uinf
	zinf 
	utau (ms)
	ztau (ms)
	gca (siemens/cm2)
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	gca = gcabar*u*u*z
	ica = gca*(v-eca)
}

INITIAL {
	rate(v)
	u = uinf
	z = zinf
}

DERIVATIVE states {
	rate(v)
	u' = (uinf-u)/utau
	z' = (zinf-z)/ztau
}
FUNCTION ualf(v(mV)) {
	LOCAL va 
	va = 19.98-v						:v+25
	if (fabs(va)<1e-04){
		va = va+0.0001
	}
	ualf = 0.19*va/(exp(va/10)-1)
}

FUNCTION ubet(v(mV)) {
	LOCAL vb 
	vb = v
	if (fabs(vb)<1e-04){
		vb = vb+0.0001
	}
	
	ubet = 0.046*(exp(-vb/20.73))
}	

FUNCTION zalf(v(mV)) {
	LOCAL va 
	va = v
	if (fabs(va)<1e-04){
		va = va+0.0001
	}
	
	zalf = 1.6e-4*exp(-va/48.4)
}

FUNCTION zbet(v(mV)) {
	LOCAL vb 
	vb = 39-v
	if (fabs(vb)<1e-04){
		vb = vb+0.0001
	}
	
	zbet = 1/(exp(vb/10)+1)
}
PROCEDURE rate(v(mV)) { LOCAL usum, zsum, ua, ub, za, zb

	ua = ualf(v) ub = ubet(v) za = zalf(v) zb = zbet(v)
	
	usum = ua+ub
	uinf = ua/usum
	utau = 1/(usum)
	
	zsum = za+zb
	zinf = za/zsum
	ztau = 1/(zsum)
}