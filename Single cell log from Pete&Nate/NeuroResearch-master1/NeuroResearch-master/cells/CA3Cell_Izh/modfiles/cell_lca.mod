:L-type voltage activated Ca current

NEURON {
    SUFFIX lca
	USEION ca READ eca WRITE ica
	RANGE gcabar, gca
	RANGE uinf, utau
}

UNITS {
    (mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcabar = 0.0001 (siemens/cm2) <0,1e9>
}

STATE { u }

ASSIGNED {
	v (mV)
	eca (mV)
	ica (mA/cm2)
	uinf
	utau (ms)
	gca (siemens/cm2)
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	gca = gcabar*u*u
	ica = gca*(v-eca)
}

INITIAL {
	rate(v)
	u = uinf
}

DERIVATIVE states {
	rate(v)
	u' = (uinf-u)/utau
}

FUNCTION ualf(v(mV)) {
	LOCAL va 
	va = 81.5-v						:v+25
	if (fabs(va)<1e-04){
		va = va+0.0001
	}
	ualf = 15.69*va/(exp(va/10)-1)
}

FUNCTION ubet(v(mV)) {
	LOCAL vb 
	vb = v
	if (fabs(vb)<1e-04){
		vb = vb+0.0001
	}
	
	ubet = 0.29*(exp(-vb/10.86))
}	

PROCEDURE rate(v(mV)) { LOCAL usum, ua, ub

	ua = ualf(v) ub = ubet(v)
	
	usum = ua+ub
	uinf = ua/usum
	utau = 1/(usum)
	
}