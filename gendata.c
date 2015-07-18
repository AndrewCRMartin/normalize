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
#include <time.h>
#include <math.h>
#include "bioplib/MathType.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXDATA 10000000
#define MAXVAL  100

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL RandomNumber(REAL maxval);


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
   int i;

   srand((unsigned int)time(NULL));
   
   for(i=0; i<MAXDATA; i++)
   {
      printf("%f String %d\n", RandomNumber((REAL)MAXVAL), i);
   }

   return(0);
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
