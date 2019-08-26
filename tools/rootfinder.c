/******************************************/
/* Root finder tools for CLASS            */
/* 20/10 2016                             */
/* Miguel Zumalacarregui                  */
/******************************************/

#include "rootfinder.h"


//TODO: add functionality for complex roots!

int rf_solve_poly_2(double a,
		   double b,
	           double c,
	           double *solutions,
		   int *complex
		   ){
  
  double x1,x2,q;
  
  if (a==0){
    
    *complex = _FALSE_;
    
    solutions[0] = -c/b;
    solutions[1] = 0;
    solutions[2] = 0;
    
    return _SUCCESS_;
    
  }

   //complex case, return real part
   if (b*b - 4.*a*c <0){
     //TODO: compiler complaints about this!
     *complex = _TRUE_;
     x1 = -b/2.;
     x2 = -b/2.;
   }
   else{
     *complex = _FALSE_;
     if (b>0)
       q = -0.5*(b + sqrt(b*b -4.*a*c));
     if (b<=0)
       q = -0.5*(b - sqrt(b*b - 4*a*c));
     
     x1 = q/a;
     x2 = c/q;
     
     if (q==0)
       x2 = 0;
  }
     
  //organize the roots smaller to larger
  if (x1<x2){
    solutions[0] = x1;
    solutions[1] = x2;
    solutions[2] = 0;
  }
  else{
    solutions[0] = x2;
    solutions[1] = x1;
    solutions[2] = 0;
  }
  
//   printf("x1 = %e, x2 = %e, complex = %i \n",solutions[0],solutions[1],complex);
  
  return _SUCCESS_;
}


/* Exact solution for real cubic polynomial
 *	sx^3 + ax^2 + bx + c = 0
 * based on the trigonometric solution
 * 
 * Fails badly when s->0, 
 * in that case use _split method
 */
int rf_solve_poly_3(double s,
		    double a,
		    double b,
	            double c,
	            double *solutions,
		    int *complex
		   ){
  
  double x0,x1,x2,q;
  double Q, R, theta, A,B;
  double x;
  int i;
  int n_sols;
  
  
  //if quadratic defer to standard case
  //TODO: consider a threshold
  if (s==0){
    rf_solve_poly_2(a,b,c,solutions,complex);
  }
  else{
    //normalize coefficients
    a /= s;
    b /= s;
    c /= s;

    // Use trigonometric formulae
    Q = (a*a -3.*b)/9.;
    R = (2.*a*a*a - 9.*a*b + 27.*c)/54.;
    
    //discriminant, case for 3 real roots
    if (R*R < Q*Q*Q){
      *complex = _FALSE_;
      
      theta = acos(R/sqrt(Q*Q*Q));
      
      x0 = -2.*sqrt(Q)*cos(theta/3.)-a/3.;
      x1 = -2.*sqrt(Q)*cos((theta+2.*_PI_)/3.)-a/3.;
      x2 = -2.*sqrt(Q)*cos((theta-2.*_PI_)/3.)-a/3.;  
      
    }//else case of complex roots
    else{
      *complex = _TRUE_;
      
      A = - pow(fabs(R) + sqrt(R*R - Q*Q*Q),1./3.);
      if (R<0)
	A *= -1.;
      
      if (A == 0)
	B = 0;
      else
	B = Q/A;
      
      x0 = (A+B) - a/3.;
      
      x1 = -0.5*(A+B)-a/3.;
      
      x2 = x1;
      
    }
    
    //TODO: organize the roots from smaller to larger
    solutions[0] = x0;
    solutions[1] = x1;
    solutions[2] = x2;
    
//     printf("x1 = %e, x2 = %e, x3 = %e, complex %i \n",x0,x1,x2,complex);
      
  }//end cubic case
  
  return _SUCCESS_;
}


/* Numerical solution for real cubic polynomial
 *	sx^3 + ax^2 + bx + c = 0
 * based on a decomposition 
 * 
 * 	sx^3 + P2(x) 
 * 
 * In the limit s->p the solutions split into
 * 
 *     x -> a/s (divergent) and x: P2(x)=0
 * 
 * The code uses those as starting points and
 * refines them numerically by Netwon's method
 * 
 * NOTE: the code assumes three solutions
 */
int rf_solve_poly_3_split(double s,
			  double a,
			  double b,
			  double c,
			  double *solutions,
			  int *complex
			  ){
  
  double f,df;
  double x0[3];
  double x;
  int i,n;
  int n_sols;
  
  //TODO: dynamical precision!
  double precision = 1e-8;
  double n_max = 100;
  
  
  //if quadratic defer to standard case
  if (s==0){
    rf_solve_poly_2(a,b,c,solutions,complex);
  }
  else{

    //TODO: separate the case with single real soluiton
    
    //put the test solutions in an array
    rf_solve_poly_2(a,b,c,x0,complex);
    // add the split solution
    x0[2] = -a/s;
    
//     printf(" x0 = %e, %e, %e \n",x0[0],x0[1],x0[2]);
    
    //run Newton's method to refine each solution
    for (i=0; i<=2; i++){
      x = x0[i];      
      f = s*x*x*x + a*x*x + b*x + c;
      n = 0;
      while (fabs(f)>precision && n < n_max){
	df = 3.*s*x*x + 2.*a*x + b;
	x -= f/df;
	//update f after df!
	f =s*x*x*x + a*x*x + b*x + c;
	n++;
      }
      solutions[i] = x;
    }
    
    //TODO: organize the roots from smaller to larger
    
//     printf("x1 = %e, x2 = %e, x3 = %e, complex %i \n",solutions[0],solutions[1],solutions[2],complex);
      
  }//end cubic case
  
  return _SUCCESS_;
}