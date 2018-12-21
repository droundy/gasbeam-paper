#include "ran.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>

#define IM1 2147483563
#define IM2 2147483563

#define AM (1.0L/IM1)
#define IMM1 (IM1-1)
#define IA1 400141L
#define IA2	406921L
#define IQ1 536681L
#define IQ2 527741L
#define IR1 122111L
#define IR2 3791L
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS DBL_EPSILON
#define RNMX (1.0L-EPS)

static long a = -1, b, c = 0, iv[NTAB];

static void array_init(void);

void init_ran(void) {
	FILE *myf = fopen("ranseed", "r");
	a = 1;
	b = 2;
	c = 3;
	if (myf) {
		if (fscanf_s(myf, " a = %ld , b = %ld , c = %ld", &a, &b, &c) == 3) {
          fclose(myf);
        } else {
          printf("Couldn't read ranseed.\n");
        }
	}
	array_init();
    return;
}
void save_ranseed() {
	FILE *myf = fopen("ranseed", "w");
	if (myf) {
		fprintf(myf, "a = %ld, b = %ld, c = %ld", a, b, c);
        fclose(myf);
	} else {
		printf("Unable to write to ranseed.  :(\n");
	}
}

void init_ran_from_file(FILE *myf) {
	a = 1;
	b = 2;
	c = 3;
	if (myf) {
		if (fscanf_s(myf, " a = %ld , b = %ld , c = %ld", &a, &b, &c) != 3) {
          printf("Error opening ran file\n");
        }
	}
	array_init();
}
void save_ranseed_to_file(FILE *myf){
	if (myf) {
		fprintf(myf, "a = %ld, b = %ld, c = %ld", a, b, c);
		
	} else {
		printf("Unable to write to file.  :(\n");
	}
}

void put_ranseed(long aa, long bb, long cc) {
	a=aa;b=bb;c=cc;
	array_init();
}
void get_ranseed(long *aa, long *bb, long *cc) {
	*aa=a;*bb=b;*cc=c;
}

static void array_init(void) {
	long j, k;
	for (j=0;j<8;j++) {
		k = a/IQ1;
		a = IA1*(a - k*IQ1) - k*IR1;
		if (a < 0) a += IM1;
		k = b/IQ2;
		b = IA2*(b - k*IQ2) - k*IR2;
		if (b < 0) b += IM2;
	}
	for (j=0;j<NTAB;j++) {
		k = a/IQ1;
		a = IA1*(a - k*IQ1) - k*IR1;
		if (a < 0) a += IM1;
		k = b/IQ2;
		b = IA2*(b - k*IQ2) - k*IR2;
		if (b < 0) b += IM2;
		iv[j] = a;
	}
}

double ran(void) {
	long j, k;
	double holder;
	if (a == -1) init_ran(); /* Someone forgot to initialize. */
	k = a/IQ1;
	a = IA1*(a - k*IQ1) - k*IR1;
	if (a < 0) a += IM1;
	k = b/IQ2;
	b= IA2*(b- k*IQ2) - k*IR2;
	if (b< 0) b+= IM2;
	j = c/NDIV;
	c = iv[j] - b;
	iv[j] = a;
	if (c < 1) c += IMM1;
	holder = AM*c;
	if (holder < 0) {
		holder = - holder;
	}
	if (holder < RNMX) {
		return holder;
	} else {
		return RNMX;
	}
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
