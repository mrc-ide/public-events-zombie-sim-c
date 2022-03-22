#ifndef MATHS_INCLUDE
#define MATHS_INCLUDE

#define SINGLE_PRECISION



#ifdef SINGLE_PRECISION
#define DTYPE float
#define DSCAN "%f\n"
#else
#define DTYPE double
#define DSCAN "%lf\n"
#endif

#include <math.h>
#include "world.h"
#include <stdio.h>
#include <stdlib.h>


class world;

#define eq(f1,f2) (fabs(f1-f2)<0.001)
#define between(x,t,b) ((x>=t) && (x<=b))
#define minF(a,b) ((a) <= (b) ? (a) : (b))
#define maxF(a,b) ((a) >= (b) ? (a) : (b))
#define PATCH_SIZE 0.166666666666667
#define CACHE_LINE_SIZE 128
#define TWO_PI 6.283185307179586476925286766559


DTYPE haversine(DTYPE lon1, DTYPE lat1, DTYPE lon2, DTYPE lat2);
DTYPE  patch_distance(world* w, int p1, int p2);
DTYPE kernel_F(world* w, double d);
int getCommunityContactPatch(world* w, double r, int patch);
double rnd(int thread_no);
long poi(double mu,int tn);
double norm(int tn);
double fsign( double num, double sign);
void initSeeds(long seed1, long seed2);
double s_expo(int tn);

#endif