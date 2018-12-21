#define PIE 3.141592653589793238462643L

double sine(double x){
	double sin;
	while(x < -1*PIE) x = x + 2*PIE;
	while(x > PIE) x = x - 2*PIE;

	if (x<=0) {
		sin = 1.27323954 * x + .405284735*x*x;
		if (sin<0) sin = .225*(sin*-sin-sin)+sin;
		else sin = .225*(sin*sin-sin)+sin;
	}
	else {
		sin = 1.27323954 * x - 0.405284735 * x * x;
		if (sin < 0) sin = .225 * (sin *-sin - sin) + sin;
		else sin = .225 * (sin * sin - sin) + sin;
	}
	return sin; 
}

double cosine(double x) {
	double sin;
	while(x < -1*PIE) x = x + 2*PIE;
	while(x > PIE) x = x - 2*PIE;

	x = PIE/2 - x;
	if (x<=0) {
		sin = 1.27323954 * x + .405284735*x*x;
		if (sin<0) sin = .225*(sin*-sin-sin)+sin;
		else sin = .225*(sin*sin-sin)+sin;
	}
	else {
		sin = 1.27323954 * x - 0.405284735 * x * x;
		if (sin < 0) sin = .225 * (sin *-sin - sin) + sin;
		else sin = .225 * (sin * sin - sin) + sin;
	}
	return sin; 
}

