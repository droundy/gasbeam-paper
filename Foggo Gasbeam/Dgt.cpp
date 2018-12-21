#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "myTrig.h"
#include <iostream>
#include "vector.h"
#include "dgt.h"

using namespace std;
void adjust_prof_for_solid_angle(double prof[90]);

vector bump_soft_sphere(vector vel, double deltatheta);
// Produces a random small vector with which to bump vel.

double pick_initial_state(vector &position, vector &velocity, 
                          double radius, double length, double nnot, double nprime, double sigma);
// Picks an initial state from off walls, from back of tube, 
// or from a hard sphere collision.  Returns the probability 
// adjustment for this state occuring.

double bin_it(double prof[90], double prob, vector position, vector velocity, double radius);
// If the particle of interest has left the tube out its front, 
// bin it in prof, and return one.  Otherwise just return zero.

void adjust_prof_for_solid_angle(double db_prof[90], double prof[90]);
// Adjusts the profile to give intensities.  Also normalizes to flux.

void bin_exp(long exp_prof[90], double costheta_for_detector, vector velocity);
// Bins the atom, taking into account the finite
// solid angle of the detector.

extern double coses[90];    // Holds the average cosine for each bin.
extern double sines[90];    // Holds the average sine for each bin.


/*==============================Main Iteratin Methods==========================================*/
void Profile(double db_prof[90], double prof[90], double radius, double length,
             double sigma, double sigmaExtra, double nnot, double nprime, double &num_counted,
             unsigned long &num_tried, long num_more,
             double costheta_for_detector, long exp_prof[90]) {
  vector position, velocity;
  // Note: Velocity should always be normalized!
  // Note: prof IS floating point.  It holds probabilities of occuring. 
  // Therefore num_counted must be floating point as well.
  double deltatheta, n, rsqr, prob, dpos;
  do{
    while ((--num_more > 0) || !num_counted) {
      dpos = radius*.5;
      ++(num_tried);
      prob = pick_initial_state(position, velocity, radius, length, nnot, nprime, 
                                sigmaExtra*.5/(1.354*1.354)/*sigma*/);
      // sigmaExtra*.5 in order to deal with CM problem with hard sphere scattering.
      // /(1.354*1.354) in order to get normal total cross sxn's back.
      rsqr = 0;
      while ( (position.z > 0) /*&& (velocity.z < 0)*/ && (rsqr < radius*radius) ) { //while it's in the tube
        dpos = radius*0.5;
        n = nnot + nprime*position.z;	//get the number density at the particles point
				
        while (nprime/n > .1/dpos) {	// Never let n drop by more than 10%.
          dpos *= 0.1;
        }
				
				
        //Soft Sphere
        deltatheta = mysqrt(n*sigmaExtra*dpos);
        //if(deltatheta > .1){
        //dpos *= .5;
        //deltatheta *= .5;
        //}
        //while(deltatheta < .06) {
        //dpos *= 2;
        //deltatheta *= 1.414214;
        //}
				
					position += velocity*dpos;
					
					if (position.z > 0) {
						velocity += bump_soft_sphere(velocity, deltatheta);
						velocity = velocity.norm();	// Normalize the velocity.
					}
					rsqr = position.x*position.x + position.y*position.y;
			}
			num_counted += bin_it(prof, prob, position, velocity, radius);
			//bin_exp(exp_prof, costheta_for_detector, velocity);
		}
	}while (num_counted < 10);
	adjust_prof_for_solid_angle(db_prof, prof);
}

void bin_exp(long exp_prof[90], double costheta_for_detector, vector velocity) {
	vector unit;
	double theta;
	unit.y = 0;
	for (theta = 0;theta < 90; ++theta) {
		unit.x = sin(theta*(PI/180));
		unit.z = -cos(theta*(PI/180));
		if (velocity*unit > costheta_for_detector) ++exp_prof[(int)theta];
	}
}

vector bump_soft_sphere(vector vel, double deltatheta) {
	// Produces a random small vector with which to bump vel.
	vector xhat(vel.z, -vel.x, vel.y), yhat; // xhat should not be parallel to vel.
	double theta = 2*PI*ran();
	xhat = cross(vel, xhat);
	xhat = xhat.norm();
	yhat = cross(xhat, vel);
	yhat = yhat.norm();
	return (xhat*cos(theta)+yhat*sin(theta))*tan(deltatheta);
}

double pick_initial_state(vector &position, vector &velocity,
				double radius, double length, double nnot, double nprime, double sigma) {
	// Picks an initial state from off walls, from back of tube, or from a hard
	// sphere collision.  Remember that the z component of velocity should be negative!
	// First decide whether to go off of a wall, a hard sphere, or in the back.
	// This preliminary stuff could be made into global stuff.
	double n1 = nnot+nprime*length,
					// The factor of two  below is due to the fact that all the ones coming in
					// the back are facing forwards.
				Pback = PI*radius*radius*n1/4, 
				Pwall = PI*radius/nprime*((n1*n1-nnot*nnot))/4 /* /2 */, 
				Phard = /*0*/sigma*PI*radius*radius/(3*nprime)*(n1*n1*n1 - nnot*nnot*nnot)/2;
	double Norm = 1/(Pwall + Phard + Pback), prob=1;
	Pwall *= Norm;
	Phard *= Norm;

	double which = ran(), r, z;
	int choice;
	vector normal(-1,0,0);
	if (which < Pwall) choice = 1;
	else if (which < Pwall + Phard) choice = 2;
	else choice = 3;
	switch(choice){
		case 1:
			// It comes off a wall.
			r = radius*.9999; 
			// Integrate: dP = ndz -> dP = n/n' dn -> P = n^2/2n'.  Then normalize by throwing
			// in factors and constants such that when P = 1, n = n1, and when P = 0, n = nnot.
			z = (sqrt((n1*n1-nnot*nnot)*ran() + nnot*nnot) - nnot)/nprime;
			velocity = ran_cos_vector(normal);	// Vel has cos dist from normal to wall.
			if (velocity.z > 0) velocity.z = -velocity.z;
			break;
		case 2:
			// It comes off a hard sphere.
			r = radius*sqrt(ran());
			// Integrate: dP = n^2dz -> dP = n^2/n' dn -> P = n^3/3n'. Then normalize by throwing
			// in factors and constants such that when P = 1, n = n1, and when P = 0, n = nnot.
			z = (pow((n1*n1*n1-nnot*nnot*nnot)*ran() + nnot*nnot*nnot, 1.0/3) - nnot)/nprime;
			velocity = ran_vector();
			if (velocity.z > 0) velocity.z = -velocity.z;
			break;
		case 3:
			// It comes in the back.
			r = radius*sqrt(ran());
			z = length;
			velocity = ran_cos_vector(normal);	// Vel has cos dist from the z axis.
		break;
	}
	position = cylindrical2vector(r, 0.0, z);
	return prob;
}

double bin_it(double prof[90], double prob, vector position, vector velocity, 
				double radius) {
	// If the particle of interest has left the tube out its front, bin it in prof, and
	// return one.  Otherwise just return zero.
	if (position.z > 0) {					// It didn't get to the end of the tube.
		return 0;
	} else if (velocity.z >= 0) {	// It is going backwards.
		return 0;
	} else {
		// We want to see if it hit the wall before leaving the tube!
		position -= velocity*(fabs(position.z/velocity.z));	
		if (position.x*position.x + position.y*position.y < radius*radius) {
			double theta = acos(-velocity.z);
			if (theta > 0) prof[(int)((180/PI)*theta)] += prob;
			if (theta < 0)prof[(int)(-(180/PI)*theta)] += prob;
			return prob;
		} else {
			return 0;
		}
	}
}

void adjust_prof_for_solid_angle(double db_prof[90], double prof[90]) {
	// Adjusts the profile to give intensities.  Also normalizes to flux.
	double flux=0;
	int i;
	for (i=0;i<90;++i) {
		flux += prof[i];
	}
	for (i=0;i<90;++i) {
		db_prof[i] = prof[i]*2*(PI/180)/(sines[i]*flux);
	}
}
