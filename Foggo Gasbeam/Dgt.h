#ifndef PI
#define PI 3.1415926
#endif

void Profile(double prof[90], double crude_prof[90], double radius, 
							double length, double sigma, double sigmaExtra, 
							double nnot, double nprime, 
							double &num_counted, unsigned long &num_tried, long num_more,
			double costheta_for_detector, long exp_prof[90]);