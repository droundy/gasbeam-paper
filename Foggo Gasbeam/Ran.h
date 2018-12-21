#include <stdio.h>

double ran();

void init_ran();
void save_ranseed();

void init_ran_from_file(FILE *);
void save_ranseed_to_file(FILE *);

void put_ranseed(long a, long b, long c);
void get_ranseed(long *a, long *b, long *c);