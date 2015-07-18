/* Error function erff() extracted from the numerics library
 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#define NUMERICS_ITMAX     100
#define NUMERICS_MAX_ERROR 5.0e-9
#define NUMERICS_FLOAT_MIN DBL_MIN
#define NUMERICS_FLOAT_MAX DBL_MAX

void NUMERICS_ERROR(const char *func, const char *msg);
double erff(double x);
double gammp(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);
void gcf(double *gammcf, double a, double x, double *gln);


void NUMERICS_ERROR(const char *func, const char *msg)
{
   fprintf(stderr, "Numerics Error in routine %s\n", func);
   fprintf(stderr, "    %s\n", msg);
}


double erff(double x)
{
   return x < 0.0 ? -gammp(.5, x * x) : gammp(.5, x * x);
}

double gammp(double a, double x)
{
   double gamser, gammcf, gln;
    
   if(x < 0.0 || a <= 0.0){
      NUMERICS_ERROR("gammp", "Invalid arguments");
      return 0.0;
   }
   
   if(x < (a+1.0)){
      gser(&gamser, a, x, &gln);
      return gamser;
   }
   else {
      gcf(&gammcf, a, x, &gln);
      return 1.0 - gammcf;
   }
}

void gser(double *gamser, double a, double x, double *gln)
{
   int n;
   double sum, del, ap;
    
   *gln = gammln(a);
   if(x <= 0.0){
      if(x < 0.0){
         NUMERICS_ERROR("gser", "x less than 0");
      }
      *gamser=0.0;
      return;
   } else {
      ap = a;
      del = sum = 1.0/a;
      for(n = 1; n <=NUMERICS_ITMAX; n++){
         ++ap;
         del *= x / ap;
         sum += del;
         if(fabs(del) < fabs(sum)*NUMERICS_MAX_ERROR){
            *gamser = sum * exp(-x+a*log(x)-(*gln));
            return;
         }
      }
   }
   NUMERICS_ERROR("gser", "a too large, NUMERICS_ITMAX too small");
   return;
}


double gammln(double xx)
{
   double x, y, tmp, ser;
   static double cof[6]={76.18009172947146, 
                         -86.50532032941677,
                         24.01409824083091, 
                         -1.231739572460166, 
                         0.1208650973866179e-2,
                         -0.5395239384953e-5};
   int j;

   y = x = xx;
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   ser = 1.000000000190015;
   for(j = 0; j <= 5; j++) ser += cof[j] / ++y;
   return -tmp + log(2.5066282746310005 * ser / x);
}

void gcf(double *gammcf, double a, double x, double *gln)
{
   int i;
   double an, b, c, d, del, h;
    
   *gln = gammln(a);
   b = x + 1.0 - a;
   c = 1.0 / NUMERICS_FLOAT_MIN;
   d = 1.0 / b;
   h = d;
   for(i = 1; i <= NUMERICS_ITMAX; i++){
      an = -i*(i-a);
      b += 2.0;
      d = an*d+b;
      if(fabs(d) < NUMERICS_FLOAT_MIN) d = NUMERICS_FLOAT_MIN;
      c = b + an / c;
      if(fabs(c) < NUMERICS_FLOAT_MIN) c = NUMERICS_FLOAT_MIN;
      d = 1.0 / d;
      del = d*c;
      h *= del;
      if(fabs(del-1.0) < NUMERICS_MAX_ERROR) break;
   }
   if(i > NUMERICS_ITMAX){
      NUMERICS_ERROR("gcf", "a too large, NUMERICS_ITMAX too small");
   }
   *gammcf = exp(-x+a*log(x)-(*gln))*h;
}
