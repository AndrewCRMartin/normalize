/*************************************************************************

   Program:    z2p
   File:       z2p.c
   
   Version:    V1.0
   Date:       24.07.09
   Function:   Calculates a 1-tailed p-value from a Z-score
   
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

   p is calculated as follows:
   ---------------------------
   D(x) = \frac{1}{2}\left[1 + 
               \mbox{erf}\left(\frac{x-\mu}{\sigma\sqrt{2}}\right)\right]
   D(x)  = \frac{1}{2}\left[1 + 
                     \mbox{erf}\left(\frac{z}{\sqrt{2}}\right)\right]
   So 
   p(|x-\mu| > |z|) = \left[1 + 
                      \mbox{erf}\left(\frac{-z}{\sqrt{2}}\right)\right]

   We then halve that to give the 1-tailed p-value

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
#define MAXBUFF 512

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL CalcProbability(REAL z);
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
   REAL zscore, p;

   if((argc != 2) || !sscanf(argv[1], "%lf", &zscore))
   {
      Usage();
      return(0);
   }
   p = CalcProbability(zscore) / (REAL)2.0;
   printf("%g\n", p);
   
   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   10.03.09 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stdout,
"\nz2p V1.0 (c) 2009, Dr. Andrew C.R. Martin, UCL\n\n\
Usage: z2p zscore\n\n\
Print's the 1-tailed p-value for a z-score to give the probability of\n\
obtaining this z-score or greater by chance. To obtain a 2-tailed\n\
p-value ensure the absolute z-score is provided and double the\n\
result.\n\n");
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



