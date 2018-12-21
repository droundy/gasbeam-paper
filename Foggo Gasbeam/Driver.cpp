#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <conio.h>
#include <iostream>
#include <fstream>
#include "dgt.h"
#include "vector.h"
#include <string.h>
#include "myTrig.h"

// #define a fun function to stop the computer.
// #ifndef __POWERPC__
// #include <conio.h>
 //#define STOP (kbhit())#if !defined(_USER32_)
// #define STOP GetKeyState(20)
// #else /* __MWERKS__ */
// #define STOP (Button())
// #endif /* __MWERKS__ */

#ifndef __POWERPC__
	#include <conio.h>
	#define STOP (_kbhit())
	//	#if !defined(_USER32_)
	//	#define STOP GetKeyState(' ')
#else /* __MWERKS__ */
	#define STOP (Button())
#endif /* __MWERKS__ */

// Must define a function for when to stop...

// Stop when the keyboard is hit...
//#define STOP GetKeyState(20)
//#define STOP (kbhit())

// Never stop!
//#define STOP (0)

using namespace std;

/* double calc_mean_vz(double crude_prof[90]); */
double calc_mean_vx(double crude_prof[90]);
double crude_prof[90], prof[90], fit[90], h1, wid1, h2, wid2;
void make_resume_file(char *fname);

double coses[90];	// Holds the average cosine for each bin.
double sines[90];	// Holds the average sine for each bin.

int main(int argc, char *argv[]) {
	int myTestCounter = 0;
	int theta;
	double nnot, nprime, sigma, radius, length, sigmaExtra,
			flow_rate, mass, temp, SS_orig, HS_orig;
	unsigned long num_tried = 0;
	double num_counted = 0;
	double costheta_for_detector;
	long exp_prof[99];

  //FILE *outFile;
  //FILE *fr_outFile;

  char outName[80], fr_name[80], yorn, junk[80];

  // First set up coses and sines to be the average value for each bin.
  for (theta=0;theta<90;theta++) {
    coses[theta] = (180/PI)*(sin((theta+1)*PI/180) - sine(theta*PI/180));
    sines[theta] = (180/PI)*(cos(theta*PI/180) - cosine((theta+1)*PI/180));
    exp_prof[theta] = 0;
  }

 // if (argc < 2) {
  cout << "What is the name of the output file? ";
  cin >> outName;

  strcpy(fr_name, outName);
  strcat(fr_name, ".fr");

  ofstream flow_rate_file;
  flow_rate_file.open(fr_name, ios::app);

  cout << "Are you continuing a previous calculation? ";
  cin >> yorn;

  /*else
    {strcpy(argv[1], outName);
    strcpy(fr_name, outName);
    strcat(fr_name, ".fr");
    //ofstream flow_rate_file(fr_name, ios::app);
    yorn = 'y';};*/

  //ofstream flow_rate_file(fr_name, ios::app);

  if (yorn == 'y') {
	//outFile = fopen(outName, "r");
   ifstream inFile(outName);
  //	if (!outFile) {
  //	  printf("Unable to open file %s\n", outName);
  //    exit(0);
  //	}

   inFile >> junk >> HS_orig;
   cout << endl << junk << ": " << HS_orig;
   inFile >> junk >> SS_orig;
   cout << endl << junk << ": " << SS_orig;
  	inFile >> junk >> sigma;
   cout << endl << junk << ": " << sigma;
	inFile >> junk >> sigmaExtra;
   cout << endl << junk << ": " << sigmaExtra;
	inFile >> junk >> radius;
	cout << endl << junk << ": " << radius;
	inFile >> junk >> length;
	cout << endl << junk << ": " << length;
	inFile >> junk >> flow_rate; // apparent flow rate
   cout << endl << junk << ": " << flow_rate;
	inFile >> junk >> flow_rate;
   cout << endl << junk << ": " << flow_rate;
   inFile >> junk >> mass;
   cout << endl << junk << ": " << mass;
	inFile >> junk >> temp;
   cout << endl << junk << ": " << temp;
	inFile >> junk >> costheta_for_detector;
   cout << endl << junk << ": " <<  costheta_for_detector;
	inFile >> junk >> nprime;
   cout << endl << junk << ": " << nprime;
	inFile >> junk >> nnot;
   cout << endl << junk << ": " << nnot;
	inFile >> junk >> num_counted;
   cout << endl << junk << ": " << num_counted;
	inFile >> junk >> num_tried;
   cout << endl << junk << ": " << num_tried;

	if(num_tried > 1L << 30) {  // Try to avoid overflow!
      cout << "Trying to avoid overflow.  Num Tried = " << num_tried << endl;
      num_tried /= 2;
      num_counted /= 2;
      cout << "Now we have Numtried = " << num_tried << " num_counted = " << num_counted << endl;
	}
	
	//fscanf(outFile, "Apparent Flow Rate\t%*lg\t%*lg\n");
	for (theta=89;theta>=0;--theta) {
      inFile >> junk >> prof[theta]  >> exp_prof[theta]  >> crude_prof[theta];
	}
//	fclose(outFile);
  } else {	// Make a new file.
	cout << "\nWhat is the hard sphere atomic cross section (Angstroms^2)? ";
	cin >> sigma;

   HS_orig = sigma;

	sigma *= 1e-8*1e-8; // convert to cross section in cm^2
	cout << "What is the extra, small angle cross section? (Ang^2 Rad^2) ";
	cin >> sigmaExtra;

   SS_orig = sigmaExtra;

	sigmaExtra *= 1e-8*1e-8;
	sigmaExtra += sigma; // I am making it all soft sphere cross section.
	sigma = 0;  //	Thus I must elimitate the hard sphere cross section.
	// Here I am increasing the sigma due to CrossSxn calculations...
	sigmaExtra *= 1.354*1.354;
	//	to get from units of Ang^2 to Ang^2 Rad^4
	cout << "What is the diameter of the tube? ";
	cin >> radius;
	radius /= 2;
	cout << "What is the length of the tube? ";
	cin >> length;
	{
      double det_dist, det_diam, sintheta;
      int i;
      cout << "What is the distance to the detector? ";
      cin >> det_dist;
      cout << "What is the diameter of the detector? ";
      cin >> det_diam;
      sintheta = (det_diam/2)/det_dist;
      costheta_for_detector = mysqrt(1 - sintheta*sintheta);
      for (i=0;i<90;i++) exp_prof[i] = 0;
	}
	cout << "What is the flow rate? ";
	cin >> flow_rate;
	cout << "What is the gas' mass (grams/mol)? ";
	cin >> mass;
    mass /= 6.0221e23;	// convert from gram/mol to grams.
	cout << "What is the temperature (Kelvin)? ";
	cin >> temp;
	double vnot = sqrt(((8.0/PI)*temp/* Kelvins */*1.3805e-16)/mass);
	// Trying new nprime and nnot stuff from Duschman
	nprime = 2*flow_rate/(vnot*radius*PI*radius*radius)*ran();
	nnot = nprime*radius*2;
  }
  double vnot = sqrt(((8.0/PI)*temp/* Kelvins */*1.3805e-16)/mass);

  // Make up a resume.bat file, before we forget.
  cout << "\nMaking resume file\n";
  make_resume_file(outName);
  // Then calculate the profile.


  /*============================================================ITERATIONS BEGIN HERE===============================================================
  =================================================================================================================================================*/
  cout << "\nIf this is a PC, press SPACE to stop.\n";

  time_t myTime = time(NULL);
  double old_fr_method, n1, d_old_fr_method;
  while (!STOP) {    
    int i;
      Profile(prof, crude_prof, radius, length, sigma,
              sigmaExtra, nnot, nprime,
              num_counted, num_tried, 370000,
              costheta_for_detector, exp_prof); //
	++myTestCounter;
	if (myTestCounter == 5) {cout << endl << "Time Ellapsed " << time(NULL) - myTime << endl;}
    // Try to avoid overflow.
    if (num_tried > 1L << 30) {  // Try to avoid overflow!
      cout << "Trying to avoid overflow.  Num Tried = " << num_tried;
      num_tried /= 2;
      num_counted /= 2;
      cout << "\nNow we have Numtried = " << num_tried << " num_counted = " << num_counted;
    }
    //////////
    n1 = nnot+nprime*length;
    old_fr_method = ((n1*n1-nnot*nnot)*PI*radius/nprime+PI*radius*radius*n1)*vnot/4*(num_counted/num_tried); //(Flux Back + Flux Wall)*(percent through)
    d_old_fr_method = old_fr_method/sqrt(num_counted);

    if ((fabs(old_fr_method - flow_rate)/d_old_fr_method > 5) &&
		(fabs(old_fr_method - flow_rate)/(old_fr_method + flow_rate) > 0.001)) {//This means 1% error	// We have a density gradient.
      double change_factor;
      cout << "\nGot " << num_counted << " through out of " << num_tried;
            cout << "\nflow rate == " << old_fr_method <<" +- " << d_old_fr_method << " nprime == " << nprime;
          //  cout << "\nprinted flow rates";
      // Remember that nnot is NOT the true nnot, but rather a measure of wall collisions.
      // Use a little viscosity to keep things under control.
     // cout << endl << "printed flow rate" << endl;
      change_factor = flow_rate/old_fr_method;	// Simply scale.
    //  cout << "\nchanged factor";
      if (change_factor > 2) {	// Someone is kookey...
        change_factor = 2;
	  } 
	  else if (change_factor < .5) {
        change_factor = 0.5;
	  }
      nprime *= change_factor;
      nnot *= change_factor;
      // Clear out all the gunk.
     // cout << "\ngunk cleared";
      num_counted = 0;
      num_tried = 0;
      cout << "\nnums reset";
      for (i=0;i<90;++i) {
        exp_prof[i] = 0;
        prof[i] = crude_prof[i] = 0;
	  }
     //cout << "\nprofile reset due to greater than 1% error";
      // Finishing setting it up to start on a new nnot.
	}
 //   cout << "\npast }";
    // Now work on the profile some more.

    // following output code moved here by DJR 8/14/00
    // this allows the program to run indefinitely (or until killed)



    //cout << "flow rate file: " << fr_name << endl;
  //  cout << endl << "about to write flow rate:" << endl;
    if ((nprime != 0) && (num_counted != 0))
	    {cout << "\nCurrent_Flow_Rate\t" << ((n1*n1-nnot*nnot)*PI*radius/(2*nprime)/2 /*4*/ + PI*radius*radius*n1/4)*vnot*num_counted/num_tried;
        flow_rate_file << endl << ((n1*n1-nnot*nnot)*PI*radius/(2*nprime)/2 /*4*/+ PI*radius*radius*n1/4)*vnot*num_counted/num_tried;}
    else cout << "\nCurrent_Flow_Rate is unavailable";


//    flow_rate_file << endl << ((n1*n1-nnot*nnot)*PI*radius/(2*nprime)/2 /*4*/+ PI*radius*radius*n1/4)*vnot*num_counted/num_tried;

   // cout << "wrote flow rate data" << endl;

    ofstream output_file(outName); //;outFile = fopen(outName, "w");
   // if (!outFile) {
   //   printf("Unable to open file %s\n", outName);
    //  exit(0);
  //  }
   // cout << "\nwriting profile to file:";
    output_file << "original_HS\t" << HS_orig;
    output_file << "\noriginal_SS\t" << SS_orig;
    output_file << "\nhard_sphere_xs\t" << sigma;
    output_file << "\nsoft_sphere_xs\t" << sigmaExtra;
    output_file << "\nradius\t" << radius;
    output_file << "\nlength\t" << length;
    //output_file << "\ncalculated_flow_rate\t" << ((n1*n1-nnot*nnot)*PI*radius/(2*nprime)/2 /*4*/ + PI*radius*radius*n1/4)*vnot*num_counted/num_tried;
  //  cout << "\nentering if statement";
    if ((nprime != 0) && (num_counted != 0))
	    output_file << "\nCurrent_Flow_Rate\t" << ((n1*n1-nnot*nnot)*PI*radius/(2*nprime)/2 /*4*/ + PI*radius*radius*n1/4)*vnot*num_counted/num_tried;
    else output_file << "\nCurrent_Flow_Rate\t" << 0;
  //  cout << "\nmade it past flow rate";
    output_file << "\nflow_rate\t" << flow_rate;
    output_file << "\ngas_mass\t" << mass;
    output_file << "\ntemperature\t" << temp;

    output_file << "\nDetector_cos_theta\t" << costheta_for_detector;

    output_file << "\nnprime\t" << nprime;
    output_file << "\nnnot\t" << nnot;

    output_file << "\nNum_Counted\t" << num_counted;
    output_file << "\nNum_Tried\t" << num_tried;

    // Then print out the profile, using both negative and positive angles.
    for (theta=89;theta>=0;--theta) {
      output_file << endl << -theta - 0.5 << "\t" << prof[theta] << "\t" << exp_prof[theta] << "\t" << crude_prof[theta];
    }
    for (theta=0;theta<90;++theta) {
      output_file << endl << theta + 0.5 << "\t" << prof[theta] << "\t" << exp_prof[theta] << "\t" << crude_prof[theta];
    }
//    fclose(outFile);
   // cout << "finished output to file" << endl;
    save_ranseed();
    //end of output code added 8/14/00 by DJR
  }
  // Someone has hit STOP.

  cout << "\nGot " << num_counted << " through out of " << num_tried;
  // We are all done.
}

double calc_mean_vx(double crude_prof[90]) {
	// I normalize within this function, so mean_vx is 1 for an orifice.
	double out=0, norm=0;
	int i;
	for (i=0;i<90;++i) {
		// The factor of 1 comes from analytically integrating this (and mean_vz) for isotropic
		// (costheta) case, and saying that they must be equal.  It comes from the splitting of
		// the velocity between vx and vy.
		out += crude_prof[i]*sines[i]/coses[i];	// tan is sin/cos.
		norm += crude_prof[i]/coses[i];
	}
	return out/(PI/4*norm);
}

void make_resume_file(char *fname) {
	ofstream resumeFile("resume.bat");
	//if (!resumeFile) {
	//	printf("Unable to create file resume.bat\n");
  //	} else {
		// Always make a backup copy!!
		resumeFile << "\ncopy " << fname << " " << fname << ".bak";
		resumeFile << "\ndgt " << fname;
		//fclose(resumeFile);
 //	}
}
