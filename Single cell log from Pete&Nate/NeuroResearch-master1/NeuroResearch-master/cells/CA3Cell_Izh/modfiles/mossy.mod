NEURON {
	POINT_PROCESS mossy
	NONSPECIFIC_CURRENT i_nmda, i_ampa
	RANGE initW
	RANGE Cdur_nmda, AlphaTmax_nmda, Beta_nmda, Erev_nmda, gbar_nmda, W_nmda, on_nmda, g_nmda
	RANGE Cdur_ampa, AlphaTmax_ampa, Beta_ampa, Erev_ampa, gbar_ampa, W_ampa, on_ampa, g_ampa
	RANGE ECa, ICa, P0, fCa, tauCa, iCatotal
	RANGE Cainf, pooldiam, z
	RANGE lambda1, lambda2, threshold1, threshold2
	RANGE fmax, fmin, Wmax, Wmin, maxChange, normW, scaleW
	RANGE pregid,postgid
	
	:Added by Ali
	RANGE U,tauRec, tauF,R,u, facfactor, fa
	RANGE aACH, bACH, aDA, bDA, wACH, wDA, calcium
}

UNITS {
	(mV) = (millivolt)
        (nA) = (nanoamp)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	pi = 3.141592 (1)
}

PARAMETER {
: parameters are vars assigned by user or changed by hoc. THey appear in nrnpointmenu
	initW = 5

	Cdur_nmda = 17.58 (ms)
	AlphaTmax_nmda = .08 (/ms)
	Beta_nmda = 0.008 (/ms)
	Erev_nmda = 0 (mV)
	gbar_nmda = .6e-3 (uS)

	Cdur_ampa = 5.31 (ms)
	AlphaTmax_ampa = 0.117 (/ms)
	Beta_ampa = 0.072 (/ms)
	Erev_ampa = 0 (mV)
	gbar_ampa = 1.7e-3 (uS)

	ECa = 120

	Cainf = 50e-6 (mM)
	pooldiam =  1.8172 (micrometer)
	z = 2

	tauCa = 50 (ms)
	P0 = .015
	fCa = .024

	lambda1 = 2.5
	lambda2 = .01
	threshold1 = 0.5 (uM)
	threshold2 = 0.6 (uM)

	fmax = 3
	fmin = .8

	:Added by Ali
	ACH = 1
	DA = 1
	LearningShutDown = 1

	tauRec = 1
	tauF = 1
	U = 0.1      :"Utilization" factor of available transmitter
	facfactor = 10
		
	aACH = 1
	bACH = 0
	wACH = 0
	aDA = 1
	bDA = 0
	wDA = 0

}

ASSIGNED {
: These are vars calculated by Neuron hoc or by the mechanism mod itself
	v (mV)

	i_nmda (nA)
	g_nmda (uS)
	on_nmda
	W_nmda

	i_ampa (nA)
	g_ampa (uS)
	on_ampa
	W_ampa

	t0 (ms)

	ICa (mA)
	Afactor	(mM/ms/nA)
	iCatotal (mA)

	dW_ampa
	Wmax
	Wmin
	maxChange
	normW
	scaleW
	
	pregid
	postgid
	
	:Added by Ali
		calcium
		u
		R
		tLast
		fa
}

STATE { r_nmda r_ampa Capoolcon }

INITIAL {
	on_nmda = 0
	r_nmda = 0
	W_nmda = initW

	on_ampa = 0
	r_ampa = 0
	W_ampa = initW

	t0 = -1

	:Wmax = fmax*initW
	:Wmin = fmin*initW
	maxChange = (Wmax-Wmin)/10
	dW_ampa = 0

	Capoolcon = Cainf
	Afactor	= 1/(z*FARADAY*4/3*pi*(pooldiam/2)^3)*(1e6)
	
	:Added by Ali
		tLast = -1e30
		u= U
		R=1-u
	fa =0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
}

DERIVATIVE release {
	if (t0>0) {
		if (t-t0 < Cdur_nmda) {
			on_nmda = 1
		} else {
			on_nmda = 0
		}
		if (t-t0 < Cdur_ampa) {
			on_ampa = 1
		} else {
			on_ampa = 0
		}
	}
	r_nmda' = AlphaTmax_nmda*on_nmda*(1-r_nmda) -Beta_nmda*r_nmda
	r_ampa' = AlphaTmax_ampa*on_ampa*(1-r_ampa) -Beta_ampa*r_ampa

	dW_ampa = eta(Capoolcon)*(lambda1*omega(Capoolcon, threshold1, threshold2)-lambda2*W_ampa)*dt

	: Limit for extreme large weight changes
	if (fabs(dW_ampa) > maxChange) {
		if (dW_ampa < 0) {
			dW_ampa = -1*maxChange
		} else {
			dW_ampa = maxChange
		}
	}

	:Normalize the weight change
	normW = (W_ampa-Wmin)/(Wmax-Wmin)
	if (dW_ampa < 0) {
		scaleW = sqrt(fabs(normW))
	} else {
		scaleW = sqrt(fabs(1.0-normW))
	}

	W_ampa = W_ampa + dW_ampa*scaleW *(1+ (wACH * (ACH - 1))) * LearningShutDown
	
	:Weight value limits
	if (W_ampa > Wmax) { 
		W_ampa = Wmax
	} else if (W_ampa < Wmin) {
 		W_ampa = Wmin
	}

	g_nmda = gbar_nmda*r_nmda
	i_nmda = W_nmda*g_nmda*(v - Erev_nmda)*sfunc(v) *R*u * facfactor :/1.5

	g_ampa = gbar_ampa*r_ampa
	i_ampa = W_ampa*g_ampa*(v - Erev_ampa)  * (aACH + (bACH * (ACH-1)))*(aDA + (bDA * (DA-1))) *R*u * facfactor :/1.5

	ICa = P0*g_nmda*(v - ECa)*sfunc(v)
	Capoolcon'= -fCa*Afactor*ICa + (Cainf-Capoolcon)/tauCa
}

NET_RECEIVE(dummy_weight) {
	t0 = t
	
	:Added by Ali, Synaptic facilitation
	:R = R*(1 - u)*exp((t-tLast)/-tauRec)+ 1 - exp((t-tLast)/-tauRec) 
	:u = u*exp((t-tLast)/-tauF)+U*(1-u*exp((t-tLast)/-tauF))
:fa = R*u * facfactor * ((R*u ) * facfactor)

	:MF
	u = 5* u * exp((t-tLast)/-tauF) + 0.05 
		
	if (u > 3) {
	u = 3
	}
	tLast = t    : Storing the most recent spike time
}

:::::::::::: FUNCTIONs and PROCEDUREs ::::::::::::

FUNCTION sfunc (v (mV)) {
	UNITSOFF
	sfunc = 1/(1+0.33*exp(-0.06*v))
	UNITSON
}

FUNCTION eta(Cani (mM)) {
	LOCAL taulearn, P1, P2, P4, Cacon
	P1 = 0.1
	P2 = P1*1e-4
	P4 = 1
	Cacon = Cani*1e3
	taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
	eta = 1/taulearn*0.001
}

FUNCTION omega(Cani (mM), threshold1 (uM), threshold2 (uM)) {
	LOCAL r, mid, Cacon
	Cacon = Cani*1e3
	r = (threshold2-threshold1)/2
	mid = (threshold1+threshold2)/2
	if (Cacon <= threshold1) { omega = 0}
	else if (Cacon >= threshold2) {	omega = 1/(1+50*exp(-50*(Cacon-threshold2)))}
	else {omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))}
}
