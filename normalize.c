/*************************************************************************

   Program:    normalize
   File:       normalize.c
   
   Version:    V1.0
   Date:       24.07.09
   Function:   Generate a normal distribution by selecting from a dataset
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2009
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      martin@biochem.ucl.ac.uk
               andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   Algorithm is as follows:
   ------------------------
   Given a set of datapoints $X$, a target mean, $\mu$, and target
   standard deviation, $\sigma$:

   foreach datapoint $X_i$
   {
      Calculate absolute $|z| = |(X_i - \mu) / \sigma|$
      Calculate probability $(p)$ of a value $(x > |z|)$
      Generate a random number $0\le r\le 1$
      if$(p \ge r)$ then place $X_i$ in the output set
   }

   p is calculated as follows:
   ---------------------------
   D(x) = \frac{1}{2}\left[1 + 
               \mbox{erf}\left(\frac{x-\mu}{\sigma\sqrt{2}}\right)\right]
   D(x)  = \frac{1}{2}\left[1 + 
                     \mbox{erf}\left(\frac{z}{\sqrt{2}}\right)\right]
   So 
   p(|x-\mu| > |z|) = \left[1 + 
                      \mbox{erf}\left(\frac{-z}{\sqrt{2}}\right)\right]

   This is taken from
      http://mathworld.wolfram.com/NormalDistribution.html
   and checked against
      http://www.stat.wvu.edu/SRS/Modules/Normal/males.html

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "erf.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXDATA 10000000
#define MAXVAL  100
#define MAXBUFF 512

typedef struct _reallist
{
   struct _reallist *next;
   REAL value;
   char data[MAXBUFF];
}  REALLIST;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REALLIST *NormalizeData(REALLIST *data, REAL targetMean, REAL targetSD);
REALLIST *ReadData(FILE *fp);
REAL CalcProbability(REAL z);
void PrintData(FILE *out, REALLIST *newdata);
REAL RandomNumber(REAL maxval);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *targetMean, REAL *targetSD);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
int main(int argc, char **argv)
{
   REALLIST *data = NULL;
   REALLIST *newdata = NULL;
   REAL targetMean;
   REAL targetSD;
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF];
   

   if(ParseCmdLine(argc, argv, InFile, OutFile, &targetMean, &targetSD))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((data = ReadData(in))==NULL)
         {
            fprintf(stderr,"Error: Unable to read input data\n");
            return(1);
         }
         if((newdata = NormalizeData(data, targetMean, targetSD))==NULL)
         {
            fprintf(stderr,"Error: Unable to build output data list\n");
            return(1);
         }

         PrintData(out, newdata);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *targetMean, REAL *targetSD)
   ----------------------------------------------------------------------
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *infile     Input filename (or blank string)
            char   *outfile    Output filename (or blank string)
            REAL   *targetMean target mean
            REAL   *targetSD   target standard deviation
   Returns: BOOL               Success

   Parse the command line

   24.07.09 Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *targetMean, REAL *targetSD)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';

   if(!argc)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2-4 arguments left                     */
         if((argc < 2) || (argc > 4))
            return(FALSE);

         /* Grab the target mean and sd                                 */
         sscanf(argv[0], "%lf", targetMean);
         argc--;
         argv++;
         sscanf(argv[0], "%lf", targetSD);
         argc--;
         argv++;

         /* Grab filenames if specified                                 */
         if(argc)
         {
            strcpy(infile, argv[0]);
         
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            
            if(argc)
               strcpy(outfile, argv[0]);
         }
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   10.03.09 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stdout,
"\nnormalize V1.0 (c) 2009, Dr. Andrew C.R. Martin, UCL\n\n\
Usage: normalize mean sd [in.dat [out.dat]]\n\n\
Samples the input dataset and writes a new set where the data are\n\
normally distributed with the required mean and standard deviation.\n");
   fprintf(stdout,
"\nThe method is to calculate the absolute Z-score of each datapoint\n\
(with respect to the required distribution) and to calculate the\n\
probability of a value having this absolute Z-score or greater.\n\
A random number between zero and one is then selected and if the\n\
probability is greater than the random number, then the datapoint\n\
is transferred to the output set.\n\n");
}


/************************************************************************/
/*>void PrintData(FILE *out, REALLIST *newdata)
   --------------------------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
void PrintData(FILE *out, REALLIST *newdata)
{
   REALLIST *n;
   for(n=newdata; n!=NULL; NEXT(n))
   {
      fprintf(out, "%s\n", n->data);
   }
}


/************************************************************************/
/*>REALLIST *NormalizeData(REALLIST *data, REAL targetMean, REAL targetSD)
   -----------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
REALLIST *NormalizeData(REALLIST *data, REAL targetMean, REAL targetSD)
{
   REALLIST *d, *n, *newdata = NULL;
   REAL z,p,r;

   srand((unsigned int)time(NULL));
   
   for(d=data; d!=NULL; NEXT(d))
   {
      z = ABS(((d->value - targetMean)/targetSD));
      p = CalcProbability(z);
      r = RandomNumber((REAL)1.0);
      if(p >= r)
      {
         if(newdata == NULL)
         {
            INIT(newdata, REALLIST);
            n = newdata;
         }
         else
         {
            ALLOCNEXT(n, REALLIST);
         }
         if(n == NULL)
         {
            FREELIST(newdata, REALLIST);
            return(NULL);
         }
         n->value = d->value;
         strcpy(n->data, d->data);
      }
   }
   return(newdata);
}


/************************************************************************/
/*>REAL RandomNumber(REAL maxval)
   ------------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
REAL RandomNumber(REAL maxval)
{
   return(maxval * rand()/(REAL)RAND_MAX);
}


/************************************************************************/
/*>REALLIST *ReadData(FILE *fp)
   ----------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
REALLIST *ReadData(FILE *fp)
{
   REALLIST *data = NULL,
      *d;
   char buffer[MAXBUFF];

   srand((unsigned int)time(NULL));
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      if(data == NULL)
      {
         INIT(data, REALLIST);
         d=data;
      }
      else
      {
         ALLOCNEXT(d, REALLIST);
      }
      if(d==NULL)
      {
         FREELIST(data, REALLIST);
      }
      sscanf(buffer, "%lf", &(d->value));
      strcpy(d->data, buffer);
   }
   return(data);
}


/************************************************************************/
/*>REAL CalcProbability(REAL z)
   ----------------------------
   Input:     
   Output:    
   Returns:   

   24.07.09  Original   By: ACRM
*/
REAL CalcProbability(REAL z)
{
   REAL p;
   p = (1+erff((double)(-z/sqrt(2))));
   return(p);
}



