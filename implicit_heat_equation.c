/*=========================================================
* 
*  NAME: implicit_heat_equation.c
*
*  This file solves the heat equation for a 1D rod using an
*  implicit method (Newton's method for iteration).
*
*  INPUTS: None
*
*  OUTPUTS: 
*
*  RETURNS:
*
*  Bryan Kaiser
*  21 August 14
*
*==========================================================*/

// Header Files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*---------------------------------------------------------*/
main(){

 // Declare variables=======================================
 int i,j,new,Nx,old;
 double beta,dx,dt,DSL,ftest,fmax,sig,time,t,T0,TL,TR;

 t = 0.0; // s, initial time
 T0 = 0.0; // C, initial rod temperature
 TL = 400.0; // C, fixed left of rod boundary temperature
 TR = 0.0; // C, fixed right of rod boundary temperature
 sig = 1.0; // m^2/s
 dx = 1.0; // m, cell length
 dt = 10; // s, time step [0.5,0.505,10]
 Nx = 50; // number of cells
 time = 100; // s, total time
 beta = 1.0/(1.0+(2.0*sig*dt)/(pow(dx, 2.0))); // denominator for Newton's method
 //printf("\nbeta = %f;\n\n",beta);
 old = 0; // timestep at n
 new = 1; // timestep at n+1 (at end of do-while loop becomes n)
 //printf("\n%d \n%d\n\n",old,new);
 ftest = 0.001*TL; // convergence criterion 
 DSL = sig*dt/(dx*dx); // diffusion stability limit

 // Computational grid=====================================

/* // Cell edges, including two ghost cells at ends
 double X[Nx+2];
 X[0] = -dx;
 for(i=1; i<(Nx+2); i++){
  X[i] = X[i-1]+dx;
 } // for i */

 // Initial conditions======================================

 // Initial temperature (w/o BCs)
 double T[2][Nx+2];
 for(i=0; i<Nx+2; i++){
  T[old][i] = T0;
  //printf("\n%f",T[old][i]); // DEBUG
 } // for i

 // Intermediate temperature array for Newton's method iteration
 double f[Nx]; // f(T) vector
 double fa[Nx]; // absolute value f(T) vector

 // Cell-Centered Computation=========================================
 
 do{
   
   // Time advancement
   t = t + dt;
   //printf("\nt=%f\n\n",t);   

   // Define BCs and update ghost zones
   // Fixes the cell edges:
   T[old][0] = 2*TL-T[old][1]; // Left BC
   T[old][Nx+1] = 2*TR-T[old][Nx]; // Right BC

   // Implicit calculation (by Newton's method)
   // 1) Initial guess for Newton's method = last step 
    for(i=0;i<Nx+2;i++){
     T[new][i] = T[old][i];
     //printf("\n%f",T[new][i]); // DEBUG
    }

    // 2) Newton's method guess for T[new][Nx]
     do {
      fmax = 0.0;
      // a) Calculate f(T):
      for(i=1;i<Nx+1;i++){
       f[i-1] = T[new][i]-T[old][i]-(sig*dt/(dx*dx))*(T[new][i+1]+T[new][i-1]-2.0*T[new][i]);  
      } // end for loop    
      // How close if f(T) to zero (i.e. calculate the 
      // absolute maximum value)?
      // b) absolute f(T)
      for(i=0;i<Nx;i++){
       fa[i] = abs(f[i]);
       //printf("\n%f",fa[i]); // DEBUG 
      } // end for loop
      // c) maximum absolute f(T)
      for(i=0;i<Nx;i++){
       if(fa[i] >= fmax){
        fmax = fa[i]; // maximum value
       } // end if 
      } // end for
      //printf("fmax = %f",fmax); // DEBUG
      // d) update T[new][Nx]
      for(i=1;i<Nx+1;i++){
        T[new][i] = T[new][i]-f[i-1]*beta;  
      } // end for loop       
     } while (fmax >= ftest); // end 2) implicit do-while loop   


   // Switch (n+1) "new" to (n) "old" for the next timestep
   for(i=1;i<Nx+1;i++){
    T[old][i] = T[new][i];
   } // end for loop 
  
 } while (t < time); // end do-while loop

 // Solution===============================================

 // Check Final Temperature Profiles
 printf("\n");
  for(i=1; i<Nx+1; i++){
    printf("\n%2.4f", T[1][i]);
  } // for i
  printf("\n\n");

 // Grid Check
 /*for(i=0; i<(Nx+2); i++){
  printf("%2.4f  ",X[i]);
 } // for i  
 printf("\n"); */ 

 printf("\nDiffusion stability condition = %2.4f\n",DSL); 
 printf("\ndx = %2.2f,dt = %2.2f,sigma = %2.2f\n",dx,dt,sig);
 printf("\nFinal time = %2.3f\n\n",t);

 // Output===================================================

 return 0;

}


