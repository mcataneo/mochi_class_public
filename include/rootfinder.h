#ifndef __RF__
#define __RF__

#define _MIN_NUMBER_OF_LAGUERRE_POINTS_ 5

/******************************************/
/* Root finder tools for CLASS            */
/* 20/10 2016                             */
/* Miguel Zumalacarregui                  */
/******************************************/
#include "common.h"

/* Structures for RF */


    /**
     * Boilerplate for C++
     */
#ifdef __cplusplus
    extern "C" {
#endif

      int rf_solve_poly_2(double a,
			  double b,
			  double c,
			  double *x,
			  int *complex
			 );
      
      //Exact trigonometric solution
      int rf_solve_poly_3(double s,
			  double a,
			  double b,
			  double c,
			  double *x,
			  int *complex
 			);
      
      //Split soluton (numerical)
      int rf_solve_poly_3_split(double s,
				double a,
				double b,
				double c,
				double *x,
				int *complex
				);
      
      

#ifdef __cplusplus
    }
#endif


#endif