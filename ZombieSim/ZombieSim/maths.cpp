#include "maths.h"

DTYPE haversine(DTYPE lon1, DTYPE lat1, DTYPE lon2, DTYPE lat2) {
  const DTYPE R = (DTYPE) 6371.297; // km
  const DTYPE PId180 = (DTYPE) (3.14159265359/180.0);
  lat1 *= PId180;
  lat2 *= PId180;
  lon1 *= PId180;
  lon2 *= PId180;
  
  DTYPE dLat = (lat2 - lat1);
  DTYPE dLon = (lon2 - lon1);
  DTYPE dLatd2 = (DTYPE) (dLat / 2.0);
  DTYPE dLond2 = (DTYPE) (dLon / 2.0);
  DTYPE a = sin(dLatd2) * sin(dLatd2) + cos(lat1) * cos(lat2) * sin(dLond2) * sin(dLond2);
  if (a > 1) a = 1; // Numerical imprecision can cause a > 1. (eg, 1 + 1E-12).
  DTYPE  c = 2 * atan2(sqrt(a), sqrt(1 - a));  // atan2 is inf when a>1, PI when a==1.
  return R * c;
}

DTYPE patch_distance(world* w, int p1, int p2) {
  if (w->patch_left_lon[p2] < w->patch_left_lon[p1]) {
    int p3 = p1;
    p1=p2;
    p2=p3;
  }

  if (eq(w->patch_left_lon[p1], w->patch_left_lon[p2])) {  // x is equal
    if (eq(w->patch_top_lat[p1], w->patch_top_lat[p2])) return 0; // Same patch
    else if (eq(w->patch_top_lat[p1], w->patch_top_lat[p2] - PATCH_SIZE)) return 0; // p2 above p1, x equal
    else if (eq(w->patch_top_lat[p1] - PATCH_SIZE, w->patch_top_lat[p2])) return 0; // p2 below p1, x equal
  } else if (eq(w->patch_left_lon[p1] + PATCH_SIZE, w->patch_left_lon[p2])) { // p1's right edge touches p2's left edge
    if (eq(w->patch_top_lat[p1], w->patch_top_lat[p2]))  return 0; // p2 on right of p1, y equal
    else if (eq(w->patch_top_lat[p1] - PATCH_SIZE,w->patch_top_lat[p2])) return 0; // p2 on right, p1 above p2 touching
    else if (eq(w->patch_top_lat[p1], w->patch_top_lat[p2] - PATCH_SIZE)) return 0; // p2 on right, p2 above p1 touching
  } 
  
  DTYPE p1_x_right = (float) (w->patch_left_lon[p1] + PATCH_SIZE);
  DTYPE p1_y_bottom = (float) (w->patch_top_lat[p1] - PATCH_SIZE);
  DTYPE p2_y_bottom = (float) (w->patch_top_lat[p2] - PATCH_SIZE);

  // Next bits: P2 always to the right of P1.
  
  // Same latitude
  if (eq(w->patch_top_lat[p1], w->patch_top_lat[p2])) return haversine(p1_x_right, p1_y_bottom, w->patch_left_lon[p2], p2_y_bottom);
  
  // P2 above (more positive) than P1
  else if (w->patch_top_lat[p2] > w->patch_top_lat[p1]) return haversine(p1_x_right, w->patch_top_lat[p1], w->patch_left_lon[p2], p2_y_bottom);
 
  // Otherwise P2 must be below (more negative than P1)
  else return haversine(p1_x_right,p1_y_bottom, w->patch_left_lon[p2], w->patch_top_lat[p2]);
 
}

DTYPE kernel_F(world* w, double d) {
  if (d > w->k_cut) return (DTYPE) 0.0;
  else return (DTYPE) (1.0 / (1.0 + (pow(d / w->k_a, w->k_b))));
}

int getCommunityContactPatch(world* w, double r, int patch) {
  unsigned int L=0;
  if (r < w->qmatrix[patch][0]) L=0;
  else if (r > w->qmatrix[patch][w->n_patches - 1]) L = w->n_patches - 1;
  else {
    unsigned int R = w->n_patches + 1;
    unsigned int P;
    double val;
    P = (R - L) / 2;
    while (P != 0) {
      P = L + P;
      if (P >= w->n_patches) {
        P = (int) w->n_patches - 1;
      }
      val = w->qmatrix[patch][P];
      if (val < r) L=P;
      else R = P;
      P = (R - L) / 2;
    }
    L++;
  }
  if (L >= w->n_patches) {
    L = (int) w->n_patches - 1;
  }
  return L;
}

// Pinched from randlib_par.cpp
    
long Xm1, Xm2, Xa1, Xa2, *Xcg1, *Xcg2, Xa1vw, Xa2vw;

long mltmod(long a,long s,long m)
/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
{
#define h 32768L
long mltmod,a0,a1,k,p,q,qh,rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
  if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
  fputs(" a, m, s out of order in mltmod - ABORT!\n",stderr);
  fprintf(stderr," a = %12ld s = %12ld m = %12ld\n",a,s,m);
  fputs(" mltmod requires: 0 < a < m; 0 < s < m\n",stderr);
  exit(1);
S10:
  if(!(a < h)) goto S20;
  a0 = a;
  p = 0;
  goto S120;
S20:
  a1 = a/h;
  a0 = a-h*a1;
  qh = m/h;
  rh = m-h*qh;
  if(!(a1 >= h)) goto S50;
  a1 -= h;
  k = s/qh;
  p = h*(s-k*qh)-k*rh;
S30:
  if(!(p < 0)) goto S40;
  p += m;
  goto S30;
S40:
  goto S60;
S50:
  p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
  if(!(a1 != 0)) goto S90;
  q = m/a1;
  k = s/q;
  p -= (k*(m-a1*q));
  if(p > 0) p -= m;
  p += (a1*(s-k*q));
S70:
  if(!(p < 0)) goto S80;
  p += m;
  goto S70;
S90:
S80:
  k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
  p = h*(p-k*qh)-k*rh;
S100:
  if(!(p < 0)) goto S110;
  p += m;
  goto S100;
S120:
S110:
  if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
  q = m/a0;
  k = s/q;
  p -= (k*(m-a0*q));
  if(p > 0) p -= m;
  p += (a0*(s-k*q));
S130:
  if(!(p < 0)) goto S140;
  p += m;
  goto S130;
S150:
S140:
  mltmod = p;
  return mltmod;
#undef h
}


double rnd(int thread_no) {
  long k, s1, s2, z;
  int curntg;
  double dev;
  curntg = CACHE_LINE_SIZE * thread_no;
  s1 = Xcg1[curntg];
  s2 = Xcg2[curntg];
  k = s1 / 53668L;
  s1 = Xa1 * (s1 - k * 53668L) - k * 12211;
  if (s1 < 0) s1 += Xm1;
  k = s2 / 52774L;
  s2 = Xa2 * (s2 - k * 52774L) - k * 3791;
  if(s2 < 0) s2 += Xm2;
  Xcg1[curntg] = s1;
  Xcg2[curntg] = s2;
  z = s1 - s2;
  if (z < 1) z += (Xm1 - 1);
  dev = ((double) z) / ((double) Xm1);
  return dev;
}

void initSeeds(long iseed1,long iseed2)
/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
  int g;
    // Max 32 threads?
  Xcg1=(long *) calloc(32*CACHE_LINE_SIZE,sizeof(long));
  Xcg2=(long *) calloc(32*CACHE_LINE_SIZE,sizeof(long));
   
  
  Xm1 = 2147483563L;
  Xm2 = 2147483399L;
  Xa1 = 40014L;
  Xa2 = 40692L;
  Xa1vw = 2082007225L;
  Xa2vw = 784306273L;
  *Xcg1 = iseed1;
  *Xcg2 = iseed2;

  for (g=1; g<32; g++) {
    *(Xcg1+(g*CACHE_LINE_SIZE)) = mltmod(Xa1vw,*(Xcg1+((g-1)*CACHE_LINE_SIZE)),Xm1);
    *(Xcg2+(g*CACHE_LINE_SIZE)) = mltmod(Xa2vw,*(Xcg2+((g-1)*CACHE_LINE_SIZE)),Xm2);
  }
}

long poi(double mu,int tn)
/*
**********************************************************************
long ignpoi_mt(double mu)
GENerate POIsson random deviate
Function
Generates a single random deviate from a Poisson
distribution with mean MU.
Arguments
mu --> The mean of the Poisson distribution from which
a random deviate is to be generated.
(mu >= 0.0)
ignpoi_mt <-- The random deviate.
Method
Renames KPOIS from TOMS as slightly modified by BWB to use RANF
instead of SUNIF.
For details see:
Ahrens, J.H. and Dieter, U.
Computer Generation of Poisson Deviates
From Modified Normal Distributions.
ACM Trans. Math. Software, 8, 2
(June 1982),163-179
**********************************************************************
**********************************************************************


P O I S S O N  DISTRIBUTION                                      


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
COMPUTER GENERATION OF POISSON DEVIATES                
FROM MODIFIED NORMAL DISTRIBUTIONS.                    
ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 

(SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  

**********************************************************************
INTEGER FUNCTION IGNPOI(IR,MU)
INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
MU=MEAN MU OF THE POISSON DISTRIBUTION
OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
SEPARATION OF CASES A AND B
*/
  {
  static double a0 = -0.5;
  static double a1 = 0.3333333;
  static double a2 = -0.2500068;
  static double a3 = 0.2000118;
  static double a4 = -0.1661269;
  static double a5 = 0.1421878;
  static double a6 = -0.1384794;
  static double a7 = 0.125006;
  /* JJV changed the initial values of MUPREV and MUOLD */
  static double fact[10] = {
    1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
    };
  /* JJV added ll to the list, for Case A */
  long ignpoi_mt,j,k,kflag,l,ll,m;
  double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
    t,u,v,x,xx,pp[35];

  if(mu < 10.0) goto S120;
  /*
  C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
  JJV changed l in Case A to ll
  */
  s = sqrt(mu);
  d = 6.0*mu*mu;
  /*
  THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
  PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
  IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
  */
  ll = (long) (mu-1.1484);
//S10:
  /*
  STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
  */
  g = mu+s*norm(tn);
  if(g < 0.0) goto S20;
  ignpoi_mt = (long) (g);
  /*
  STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
  */
  if(ignpoi_mt >= ll) return ignpoi_mt;
  /*
  STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
  */
  fk = (double)ignpoi_mt;
  difmuk = mu-fk;
  u = rnd(tn);
  if(d*u >= difmuk*difmuk*difmuk) return ignpoi_mt;
S20:
  /*
  STEP P. PREPARATIONS FOR STEPS Q AND H.
  (RECALCULATIONS OF PARAMETERS IF NECESSARY)
  .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
  THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
  APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
  C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
  */
  omega = 0.3989423/s;
  b1 = 4.166667E-2/mu;
  b2 = 0.3*b1*b1;
  c3 = 0.1428571*b1*b2;
  c2 = b2-15.0*c3;
  c1 = b1-6.0*b2+45.0*c3;
  c0 = 1.0-b1+3.0*b2-15.0*c3;
  c = 0.1069/mu;
//S30:
  if(g < 0.0) goto S50;
  /*
  'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
  */
  kflag = 0;
  goto S70;
S40:
  /*
  STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
  */
  if(fy-u*fy <= py*exp(px-fx)) return ignpoi_mt;
S50:
  /*
  STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
  DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
  (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
  */
  e = s_expo(tn);


  u = rnd(tn);
  u += (u-1.0);
  t = 1.8+fsign(e,u);
  if(t <= -0.6744) goto S50;
  ignpoi_mt = (long) (mu+s*t);
  fk = (double)ignpoi_mt;
  difmuk = mu-fk;
  /*
  'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
  */
  kflag = 1;
  goto S70;
S60:
  /*
  STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
  */
  if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
  return ignpoi_mt;
S70:
  /*
  STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
  CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
  */
  if(ignpoi_mt >= 10) goto S80;
  px = -mu;
  py = pow(mu,(double)ignpoi_mt)/ *(fact+ignpoi_mt);
  goto S110;
S80:
  /*
  CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
  A0-A7 FOR ACCURACY WHEN ADVISABLE
  .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
  */
  del = 8.333333E-2/fk;
  del -= (4.8*del*del*del);
  v = difmuk/fk;
  if(fabs(v) <= 0.25) goto S90;
  px = fk*log(1.0+v)-difmuk-del;
  goto S100;
S90:
  px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
S100:
  py = 0.3989423/sqrt(fk);
S110:
  x = (0.5-difmuk)/s;
  xx = x*x;
  fx = -0.5*xx;
  fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
  if(kflag <= 0) goto S40;
  goto S60;
S120:
  /*
  C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
  JJV changed MUPREV assignment to initial value
  */
  m = maxF(1L,(long) (mu));
  l = 0;
  p = exp(-mu);
  q = p0 = p;
S130:
  /*
  STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
  */
  u = rnd(tn);
  ignpoi_mt = 0;
  if(u <= p0) return ignpoi_mt;
  /*
  STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
  PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
  (0.458=PP(9) FOR MU=10)
  */
  if(l == 0) goto S150;
  j = 1;
  if(u > 0.458) j = minF(l,m);
  for(k=j; k<=l; k++) {
    if(u <= *(pp+k-1)) goto S180;
    }
  if(l == 35) goto S130;
S150:
/*
STEP C. CREATION OF NEW POISSON PROBABILITIES P
AND THEIR CUMULATIVES Q=PP(K)
*/
  l += 1;
  for(k=l; k<=35; k++) {
    p = p*mu/(double)k;
    q += p;
    *(pp+k-1) = q;
    if(u <= q) goto S170;
    }
  l = 35;
  goto S130;
  S170:
  l = k;
  S180:
  ignpoi_mt = k;
  return ignpoi_mt;
}

double norm(int tn)
/*
**********************************************************************


(STANDARD-)  N O R M A L  DISTRIBUTION                           


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
SAMPLING FROM THE NORMAL DISTRIBUTION.                 
MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
(M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
SUNIF.  The argument IR thus goes away.                          

**********************************************************************
THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
  {
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
      0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
      0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
      1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
      1.862732,2.153875
    };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
      0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
      0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
      0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
    };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
      1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
      2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
      4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
      9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
    };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
      4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
      4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
      5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
      8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
    };
  long i;
  double snorm_mt,u,s,ustar,aa,w,y,tt;
  
  u = rnd(tn);
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
  START CENTER
  */
  ustar = u-(double)i;
  aa = *(a+i-1);
S40:
  if(ustar <= *(t+i-1)) goto S60;
  w = (ustar-*(t+i-1))**(h+i-1);
S50:
  /*
  EXIT   (BOTH CASES)
  */
  y = aa+w;
  snorm_mt = y;
  if(s == 1.0) snorm_mt = -y;
  return snorm_mt;
S60:
  /*
  CENTER CONTINUED
  */
  u = rnd(tn);
  w = u*(*(a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = rnd(tn);
S80:
  if(ustar > tt) goto S50;
  u = rnd(tn);
  if(ustar >= u) goto S70;
  ustar = rnd(tn);
  goto S40;
S100:
  /*
  START TAIL
  */
  i = 6;
  aa = *(a+31);
  goto S120;
S110:
  aa += *(d+i-1);
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u**(d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = rnd(tn);
  if(ustar > tt) goto S50;
  u = rnd(tn);
  if(ustar >= u) goto S150;
  u = rnd(tn);
  goto S140;
}

double fsign( double num, double sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
else return num;
}

double s_expo(int tn)
/*
**********************************************************************


(STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                


**********************************************************************
**********************************************************************

FOR DETAILS SEE:                                                 

AHRENS, J.H. AND DIETER, U.                            
COMPUTER METHODS FOR SAMPLING FROM THE                 
EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               

ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       

Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
SUNIF.  The argument IR thus goes away.                          

**********************************************************************
Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
(HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
  {
  static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
      .9999999
    };
  long i;
  double sexpo_mt,a,u,ustar,umin;

  a = 0.0;
  u = rnd(tn);
  goto S30;
S20:
  a += q[0];
S30:
  u += u;
  /*
  * JJV changed the following to reflect the true algorithm and prevent
  * JJV unpredictable behavior if U is initially 0.5.
  *  if(u <= 1.0) goto S20;
  */
  if(u < 1.0) goto S20;
  u -= 1.0;
  if(u > q[0]) goto S60;
  sexpo_mt = a+u;

  
  return sexpo_mt;
S60:
  i = 1;
  ustar = rnd(tn);
  umin = ustar;
S70:
  ustar = rnd(tn);
  if(ustar < umin) umin = ustar;
  i += 1;
  if(u > q[i-1]) goto S70;
  
  
  return  a+umin*q[0];
}
