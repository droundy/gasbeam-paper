#include <math.h>
#include "myTrig.h"


#ifndef PI
#define PI 3.141592653590L
#endif

// Random number generator.  Requires ran.c for use.

double ran(void);			// Automatically initializes the first time, and saves ranseed at exit.
void save_ranseed();	// Manually save ranseed.

// Utility functions.

inline double myatan2(double x, double y);	// A non-crashable atan2. (y may be zero)
inline double mysqrt(double x);							// A non-crashable sqrt. (x may be negative)

// A fast approximate sin and cos.
//void mytrig(double theta, double *costheta, double *sintheta);

double myexp(double x);

// Vector functions.

struct vector {
	inline vector(void);
	inline vector(double xx, double yy, double zz);
	double x, y, z;
	inline double operator*(vector b);
	inline vector operator*(double k);
	inline vector operator/(double k);
	inline vector operator+(vector b);
	inline vector operator-(vector b);
	inline vector operator*=(double k);
	inline vector operator/=(double k);
	inline vector operator+=(vector b);
	inline vector operator-=(vector b);
	inline double mag(void);
	inline vector norm(void);
};

inline vector cross(vector a, vector b);
inline double dot(vector a, vector b);
inline vector polar2vector(double r, double theta, double phi);
inline vector cylindrical2vector(double r, double theta, double z);
inline void vector2polar(vector v, double *r, double *theta, double *phi);
inline void vector2cylindrical(vector v, double *r, double *theta, double *z);
inline vector ran_vector(void);	// Picks a randomly oriented unit vector.
inline vector ran_cos_vector(vector direction);	// Picks a randomly oriented unit vector.

////////////////////// The Implementations are below! //////////////////////////////////

inline double myatan2(double x, double y) {
	if (y) {
		return atan2(x,y);
	} else {
		return PI/2;
	}
}

inline double mysqrt(double x) {
	if (x>=0) {
		return sqrt(x);
	} else {
		return sqrt(-x);
	}
}

inline vector::vector(void) {
	x=0;y=0;z=0;
}

inline vector::vector(double xx, double yy, double zz) {
	x=xx;y=yy;z=zz;
}

inline double vector::operator*(vector b) {	// Dot product
	return x*b.x + y*b.y + z*b.z;
}
	
inline vector vector::operator*(double k) {
	return vector (k*x,k*y,k*z);
}
	
inline vector vector::operator/(double k) {
	return vector (x/k,y/k,z/k);

}
	
inline vector vector::operator+(vector b) {
	return vector (b.x+x,b.y+y,b.z+z);
}
	
inline vector vector::operator-(vector b) {
	return vector (x-b.x,y-b.y,z-b.z);
}
	
	
inline vector vector::operator*=(double k) {
	x *= k; y *= k; z *= k;
	return *this;
}
	
inline vector vector::operator/=(double k) {
	x /= k; y /= k; z /= k;
	return *this;
}
	
inline vector vector::operator+=(vector b) {
	x += b.x; y += b.y; z += b.z;
	return *this;
}
	
inline vector vector::operator-=(vector b) {
	x -= b.x; y -= b.y; z -= b.z;
	return *this;
}
	
inline double vector::mag(void) {
	return sqrt(x*x+y*y+z*z);
}
	
inline double abs(vector v) {
	return v.mag();
}

inline vector vector::norm(void) {
	vector out(x, y, z);
	return out / out.mag();
}

inline vector cross(vector a, vector b) {
	return vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline double dot(vector a, vector b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline vector polar2vector(double r, double theta, double phi) {
	return vector (r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
}

inline vector cylindrical2vector(double r, double theta, double z) {
	return vector (r*cos(theta), r*sin(theta), z);
}

inline void vector2polar(vector v, double *r, double *theta, double *phi) {
	*r = v.mag();
	*theta = acos(v.z/ *r);
	*phi = myatan2(v.x, v.y);
}

inline void vector2cylindrical(vector v, double *r, double *theta, double *z) {
	*r = sqrt(v.x*v.x+v.y*v.y);
	*theta = myatan2(v.x, v.y);
	*z = v.z;
}

inline vector ran_vector(void) {
	double costheta = 2*ran()-1, sintheta = mysqrt(1-costheta*costheta), phi=(2*PI)*ran();
	return vector (sintheta*cos(phi), sintheta*sin(phi), costheta);
}

inline vector ran_cos_vector(vector direction) {
	double cossqrtheta = ran(), sintheta = mysqrt(1 - cossqrtheta), costheta = mysqrt(cossqrtheta);
	double phi = 2*PI*ran();
	vector xhat(direction.z, -direction.x, direction.y), yhat, out;	// Pick a basis set.
	xhat = cross(xhat, direction);
	xhat = xhat.norm();
	yhat = cross(direction, xhat);
	return direction*costheta + xhat*(sintheta*cos(phi)) + yhat*(sintheta*sin(phi));
}