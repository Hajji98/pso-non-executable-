/***************************************************************************
Some binary PSO derivations
  First version          : 2004-10
    Last update          : 2005-02-02
    email                : Maurice.Clerc@WriteMe.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

 /*


  Special Particle Swarm Optimiser for binary problems.

  Actually, this program contains several algorithms (see the "method" parameter)
  founded on the PSO canonical model.
  The purpose of this program is indeed to show how easy it is to derive some good
  algorithms from this model.

  The "method" parameter is an arbitrary code.
  If method<100
	information links are modified at random after each iteration
   if there has been no improvement. The swarm size is S.
   Each particle choose at random(K choices) which ones it informs. 
   It means each particle is informed by a number of particles
   between 0 and S, but, of course, probably more or less equal to K

  If method>100
   Methods originally designed with information links  that use
   the classical circular topology.
   Each particle has then always exactly the same number K of informers
   (including itself).


  You may use a pseudo-gradient option to find the best informer
  (using a cyclic distance between binary positions).
  Note that is does really make sense only for some methods (like 1)

  Some binary methods are based on two cyclic algebras:
   - in {-1,0,1} for the velocity
   - in {0,1} for the position
  The similarity with a cellular automaton is then particularly obvious
  when you use just one particle.

  You will see that some methods like 0 are _too_ good on problems 1,2 and 3.
  Although they are often called "deceptive" problems, they are in fact,
  if I dare say, just "deceptive deceptive" ones:
  the structure of the solution is so particular that they are indeed extremely easy.
  They can be solved very quickly with a "swarm" of size one!

  So, I have added some other problems, in particular Zebra3.

  The source code is a bit complicated, for I have added the possibility
  to combine several method in an adaptive way (parameter adapt_method).

  My advice is to set it first to 0, and to try individually different methods.
  Unlike in TRIBES, S and K are not adaptive (see TO DO, though)
  so you have to choose them manually. 

  I have added a "quality measure" that gives an estimation of how good is a given
  algorithm on a given problem, assuming you run it in a loop
  with different maximum evaluations
  (see eval_max1, eval_max2, delta_eval_max)

  For some functions, I have also added a trick to save some computational time:
 when the function value is progressively computed as a sum of non negative values
 there is usually no need to compute it entirely. You just have to be sure it does
 not improve the best solution known by the particle

  Have fun and keep me posted.

  P.S. 1. If your compiler uses a reasonably good pseudo random numbers generator
     you don't have to use KISS.

  P.S. 2. Some options and subprogramms are not really useful. 
       Don't forget this is a research version!

  P.S. 3 There is an almost separate subprogram (fundamental_hyp)
          to perform some statistical tests about the hypothesis "nearer is better"

 */

/* TO DO
-	In fact, there _is_ an adaptive K option, but far too rudimentary (either 2 or 3).
	Has to be improved.
- add Parisian method for some kinds of problems?

*/

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <time.h>

#define	D_max 150  // Max number of dimensions of the search space
#define E_max 1000 // Max number of executions for one problem
#define	S_max 200 // Max swarm size
#define F_max 3 // (for multiobjective) Max number of functions
#define ulong unsigned long 	// For KISS pseudo random numbers generator
#define RAND_MAX_KISS ((unsigned long) 4294967295) // idem
					// You don't need KISS if you have a good rand() function
#define MAX_PEAKS 500

// Structures
struct f	{int size;double fi[F_max];}; // Fitness's
struct	bits	{int size;int b[D_max];};
struct	position {struct bits x; struct f f;};


// Sub-programs
double	alea(double a,double b);
int		alea_integer(int a,int b);
double alea_normal(double mean, double std_dev);
struct bits binary(unsigned int a,int size);
int		best_(struct f fx,struct f fy);
int		best_info(int s,int option);
double combin(int D,int k);
int diameter();   // Swarm diameter
double	distance(int s1,int s2, int option);
void  fundamental_hyp(int function);
void init_swarm(int option);
double number(int d,struct bits x,double scale);
double	perf(int s,int function);
ulong	rand_kiss();
void	seed_rand_kiss(ulong seed);
void swarm_size(); // Just for plotting some curves
double	total_error(struct f f);

// Global variables
struct	position best;	// Best of the best position
int	D;				// Search space dimension
int Sp; // pseudo swarm size for problem 11
int Diam; // Swarm diameter
double  E; // Useful for some test functions
struct f	fMin;	// Objective(s) to reach
double	granul[D_max];	// Granularity for each dimension (0 = continue, 1 = integer)
int init; // Flag. Just to know wether evaluations are in init phase or not
int		nb_eval; // Total number of evaluations
int	n_exec,n_exec_max; // Nbs of executions
int	LINKS[S_max][S_max]; // Information links
double mean_val;

int     n_iter;
struct	position P[S_max];	// Positions
struct	position P_m[S_max];// Best positions found by each particle
struct position P_init[S_max];
int	S;				// Swarm size
double pi; // Useful for some test functions
int pivots[S_max]; // Flags to know which particles are pivots
double  Q[D_max][D_max]; // Matrix for quadratic problems
double save_time; // To evaluate computational time saved
                // by not entirely computing the fitness
struct	bits V[S_max];// Velocities
double	xmin[D_max],xmax[D_max];	// Intervals defining the search space

// For Multimodal random generator
int peaks;                          /* The number of peaks */
int landscape[MAX_PEAKS][D_max]; /* Holds the location of the peaks */

// File(s)
FILE	*f_multi;	// In case of multiobjective, save results, in order to later find
					// the Pareto front
FILE *f_Q; // Matrix for quadratic problems
FILE *f_Qr; // Matrix for random quadratic problems
FILE *f_landscape; // Landscape for Multimodal problem
FILE *f_init; // If the initial position is read from a file
FILE *f_trace;
FILE *f_SSS; // Structured search space estimation. Just for information
			// If the parameter SSS is set to 1, this file is generated
			// and there is no optimisation
FILE *f_synth; // Summary of the results

void main()
{
int	adapt_K;
double adapt_latency;
int	adapt_method;
double alpha;
struct position best_prev, best_temp;
int bit,bit0,bit1,bit2,bit3;
struct bits B_mod[S_max]; // Modified bits
double	c2,c3;// Confidence coefficients
double cp,cp0,cp1,cp2,cp3;
double c_pr0;
int check_hyp;
int choice[3];
int	d;				// Current dimension
int DD;
double density;
int deter;  // Flag for deterministic or random search
int distrib[D_max]; // Just for some statistics
int dist;
double	eps;		// Wanted accuracy
double	eps_mean;	// Average accuracy
double	error;		// Total error for a given position
double	error0,error1,error2, error3;
double	error_prev; // Total error after previous iteration
double error_prev2;
double error_[E_max];
int	eval_max;		// Max number of evaluations
int eval_max1, eval_max2, delta_eval_max;
double	eval_mean;	// Mean number of evaluations
int fixed_link; 
int	function[F_max];
int	g;				// Rank of the best informer
//int g_[S_max]; // Best informers list
struct position Gp;
int	i,j;
int	init_links;		// Flag to (re)init or not the information links
int init_option;  // 0 => as usually, at random  for the positions
                  // > 0. See  init_swarm(), 
                  // some attempts to find a more evenly distribution
                  // <0 => read on the file init.txt
int	k,k0, k1,k2;
double k_tot[D_max]={0}; // Just for some distribution tests
int	K;				// Max number of informed particles by a given one
double K_choice[3][2];
//int	K_prev;
int	m;
int method;
double min; // Best result along several runs
int n1,n2,n3;
int n_adapt,n_adapt_max;
int n_adaptK;
int nb_pivots;
int		n_f;		// (multiobjective) number of functions to optimize
int	n_failure;		// Number of failures

int	 n_iter0,n_iter1; // Number of iterations
int method_rank;
int	option_best; // Option for the best informer choice
						// 0 => really the best (classical)
						// 1 => using pseudo-gradient
struct position P_prev[S_max], P_temp[S_max];
struct position P_m_prev[S_max], P_m_temp[S_max];
struct position Pp;
double  pr; 
double pr0; // Strategy probabilities           
int print_level; // > 0 => more verbose
struct bits V_prev[S_max], V_temp[S_max];
struct position Vp;
double quality, quality_ref;
double r;
int	s,s1,s2; // Rank of the current particle
int SSS; // Structured Search Space option
int strategies[5];
int strat_number;
double success_rate;
double t1,t2;
int tribe; // 0 => classical. after the move, update just the previous best
          // of the particle
            // #0 => update the whole tribe
//float z;
double z;
float zz;
double z7;


f_synth=fopen("f_synth.txt","w");
f_trace=fopen("f_trace.txt","w");


//swarm_size(); goto end; // Just to plot some curves

seed_rand_kiss(1); // Initialise KISS
E=exp(1);
pi=acos(-1);

// ------------------------------------------------------ Control parameters
option_best=0;	// 1 or 0. Using pseudo-gradient or not
init_option=2; // 0 => positions completely at random
              // 1  (removed. Was not good)
              // 2 => semi-deterministic distribution
              
			// <0 => read from file init.txt
tribe=0; // 0 => classical. After a move, update just the previous best
        // of the current particle
        // 1 => update for the whole current tribe
        // 2 => update the worst of the current tribe
        // WARNING: here "tribe" means "informant group"
				
method=11; // cf MOVES
		 // Method 0 is almost a hoax. Try to guess why ;-)
   // Method 6 is simple and _sometimes_  good 
   // Method 7 seem to be globally the best (or method 16?)
   // Have a look at method 8, particularly for quadratic problems
	// Method 11 is quick on easy problems like 1, 2, 4 etc.
  // Method 16 is like 11 but with "taboo". Usually works better
	// Method 100 is for comparison (Kennedy's method, improved by me)


S=12;	// Swarm size. Typically 35, but just 2 may be enough for particular problems
      // like 1, 2, 3, 
      // On the other hand, you may need more particles for high dimension problems (>40)
      // For the "optimal" swarm size, have a look at the paper on my PSO site:
      // a rough approximation of the theoretical formula is 
		//  S=INT(9.5+ 0.124*(Dimension-9))

K=3;	// Max number of particles informed by a given one
      // Usually 2 or 3 for difficult problems
      // although K=1 sometimes works 
	 // if <0, K=1 or 3
	// For method >= 100, you may use K=S
fixed_link=method>=100;   // Mainly for methods>=100. But you may try with some others

adapt_method=0;
	// 0 => always the same method
	// 1 => sequential mode, if there has been no improvement for a while, tries another method
	// 2 => parallell mode: from time to time tries TWO methods and keep the best

if (adapt_method>0) // Strategy set definition. The previous "method" value is ignored
{
	strat_number=2; 
	strategies[0]=16;
	strategies[1]=0;

	strategies[2]=5;
	strategies[3]=5;
	strategies[4]=5;
}

deter=0; // 0 => some randomness, 1 => deterministic (for binary methods)
        // -1 => more randomness
        // WARNING: "old" parameter. Not valid for all methods
        
print_level=-1; // Just more or less verbous
//--------------------------------------------------------- Problem
D=30;	// Search space dimension
n_f=1;			// Nb of functions to simultaneously optimise
function[0]=1; //Function code (see perf())
Sp= 10; // Pseudo swarm size for problem 11
SSS=0; // number of sampling points. If >0 => just write the SSS file. No optimisation

/*
1 Goldberg's order 3 deceptive problem (binary problem)
2 Bipolar order 6
3 Muhlenbein's order-5
4 Clerc's order 3 1
5 Clerc's order 3 2
6 Clerc's order 3 (Zebra3)
7 Parabola
8 Using Multimodal bitstring random generator. Use 100D, 20 peaks for test
9 Clerc's Trap
10 Whitney DF2
11 Look for a good initialisation for further problems (see perf() )
-1 Quadratic (data on a file)
-2 Quadratic (random matrix)
-8 Multimodal (data on a file)

99 Test 
*/


//--------------------------------- 
eps=0.00; 		// Accuracy. Can be null, for it is a discrete problem
fMin.size=n_f;
for (i=0;i<n_f;i++)
	fMin.fi[i]=0.0;	// Objective(s) to reach
n_exec_max=100; 		// Numbers of runs
eval_max1=20000; 	// Min Max number of evaluations
eval_max2=20000;   // Max Max number of evalutions
delta_eval_max=1000;  // Step
  density=1.0; // For Quadratic problem -2
  peaks=20; // For Multimodal 8

check_hyp=0; //1 =>Separate test, for the fundamental hypothesis
            // If set to 1, NO optimisation at all
            
//========================================================= Algorithm
if(K<0) adapt_K=1;else adapt_K=0; K=abs(K);
n_adapt_max=(int)((double)S/(double)(K-1)); // Number of evaluations between two adaptations

if (function[0]==-1)   // Mono objectif xQx problem. Q on a file

{
  	f_Q=fopen("Q100paper.txt","r"); // Read the matrix
  // 	f_Q=fopen("Q12.txt","r"); // Read the matrix  
   printf("\n Q read from file");
   fscanf(f_Q,"%i",&D);  // WARNING. Read the dimension on the file
   printf("\n Q matrix %ix%i",D,D);


  density=0;
   for (i=0;i<D;i++)
    for (d=0;d<D;d++)
      {
        fscanf(f_Q,"%f",&zz);
       Q[i][d]=zz;
       if(fabs(zz)>0) density=density+1;
       }
       z=D;
    density=density/(D*D);
    printf("  Density %.2f",density);
    printf("\n Minimum value 0");
    printf("\n");
}
if(function[0]==-2) // Mono objectif xQx problem. Random Q
{
  printf("\n Random Q matrix %ix%i",D,D);
  printf("  Density %f",density);
  printf("\n Unknown minimum value.  ...");

  printf("\n");
  f_Qr=fopen("Qr.txt","w"); // Save the matrix
  printf("\n Q saved on Qr.txt");
  fprintf(f_Qr,"%i",D);
  for (i=0;i<D;i++)
  {
    fprintf(f_Qr,"\n");
    for (d=0;d<D;d++)
    {
       if(alea(0,1)<density)
       n1= alea_integer(-3*D,3*D);
        Q[i][d]=n1;
        fprintf(f_Qr,"%i ",n1);
    }
  }
  function[0]=-1;
}


if (function[0]==8) // Multimodal bitstring random generator
                    // (William Spears)
{ 
f_landscape=fopen("f_landscape.txt","w");

 for(i = 0; i < peaks; i++)
  {
       for(j = 0; j < D; j++)
        {
	        landscape[i][j] = rand()&01;
       }
  }

fprintf(f_landscape,"%i %i",D,peaks);  
for (i=0;i<peaks;i++)
{
  fprintf(f_landscape,"\n");

  for (j=0;j<D;j++)
  fprintf(f_landscape,"%i ",  landscape[i][j]);
}
}

if (function[0]==-8) // Multimodal bitstring random generator
{ 
  f_landscape=fopen("f_landscape_12_20.txt","r");
  fscanf(f_landscape,"%i %i", &D, &peaks);
  for(i = 0; i < peaks; i++)
  {
       for(j = 0; j < D; j++)
        {
	        fscanf(f_landscape,"%i",&landscape[i][j]);
       }
  }
  function[0]=8;
 }

 
if (SSS==0) goto check_h;
printf("\n Structured Search Space on file SSS.txt");
f_SSS=fopen("SSS.txt","w");
fprintf(f_SSS,"Distance_to_solution  Fitness");

printf("\nDistance_to_solution  Fitness");
P[0].x.size=D;

switch (function[0])
{
default:
for (s=0;s<D;s++) P[0].x.b[s]=0;
break;


case 1:
for (s=0;s<D;s++) P[0].x.b[s]=1;
break;


case 6:

m=0;
for (d=0;d<D;d=d+3)
{
  m=m+1;
  if(m%2==0)

  {
    for (i=d;i<d+3;i++) P[0].x.b[d]=1;
   }
   else
   {
    for (i=d;i<d+3;i++) P[0].x.b[d]=0;
   }

}

break;
}

P[1].x.size=D;
P[1].f.size=n_f;

for (m=0;m<SSS;m++)
{

		for (d=0;d<D;d++) P[1].x.b[d]=alea_integer(0,1);

		for (i=0;i<n_f;i++)
			P[1].f.fi[i]=fabs(perf(1,function[i])-fMin.fi[i]);
		error=total_error(P[1].f);
	// Compute distance to solution
	r=distance(0,1,-1);	
printf("\n %i %.0f %f",m,r,error);
fprintf(f_SSS,"\n%.0f %f",r,error);
}

goto end;

check_h:
if(check_hyp==0) goto optim;
  fundamental_hyp(function[0]);
  goto end;

//----------------------------------------------- OPTIMISATION

optim:
t1=clock();
mean_val=pow(2.0,D-1);

// Save if multiobjective
if (n_f>1)
	f_multi=fopen("f_multi.txt","w");

alpha=pow(0.99,1.0/(3*D-1)); // This is just to estimate the "quality" of the algorithm.
							// Better to run on several eval_max.
							// For the "ideal" algo, we want a succes rate of 0.99
							// after at most 3 evaluations for each dimension

quality_ref=0;

//------------------------------LOOP on maximum number of evaluations
for (eval_max=eval_max1;eval_max<=eval_max2 ;eval_max=eval_max+delta_eval_max) 
{

// Initialisation of information variables
n_exec=0; eval_mean=0; eps_mean=0; n_failure=0;
save_time=0;

// x_p = best previous position
// x_g = best previous position of the best informer
cp1=1.;
if(deter==1) {cp2=1; cp3=1;} // Deterministic: always towards x_p and x_g

init:          // Initialisations of positions and velocities
init=1;
n_iter=-1;
n2=1; n3=1; // For method 7

if (adapt_method>0) // If adaptive method
{
	n_adapt=0; method_rank=0;
	method=strategies[0];
	adapt_latency=1;
	n_iter0=n_adapt_max;
	n_iter1=0;
	n_iter1=n_adapt_max;
// if(method==100) K=35; else K=alea_integer(2,3);   //TEST
}


if (adapt_K==1) // If adaptive K
{
	K=1+2*alea_integer(0,1); // 1 or 3
	K_choice[0][0]=0; K_choice[2][0]=0; // Cumulated improvement
	K_choice[0][1]=0; K_choice[2][1]=0; // Number of times K is used
	n_adaptK=0;
}

// Some confidence coefficients for some methods
pr0=0.5+0.5/log((1+log((double)D))); // a
//pr0=1/pow(2,D); // b
c_pr0=0.5*log((double)D)/log((double)S);

n_exec=n_exec+1;


  // Initialise particles and velocities
  init_swarm(init_option);

// First evaluations
nb_eval=0;

for (s=0;s<S;s++)
{
  P[s].f.size=n_f;
	for (i=0;i<n_f;i++)
		P[s].f.fi[i]=fabs(perf(s,function[i])-fMin.fi[i]);
	P_m[s]=P[s]; // Best position = current one
}



// Save the best of the bests
best=P_m[0];
for (s=0;s<S;s++) if (best_(P_m[s].f,best.f)==1) best=P_m[s];


// Current min error
error=total_error(best.f);
error_prev=error;
error_prev2=error_prev;
init_links=1; // So that information links will be initialized

// For method 7
cp0=1-0.5*log((double)D)/log((double)S);
if (K>1) cp0=cp0*log((double)K);
 //z7=(double)nb_eval*K/(double)S;
  z7=(double)K/(double)S;
 if (z7>1) z7=log(z7);

 cp1=cp0/z7;
 cp2=cp1;

 if(K>2)
  cp3=cp1/log(K-1);// Decreasing. 
 else cp3=cp1;

   n2=(int)(cp2*D);// Decreasing. The more you trust p, the less you update it
   n3=(int)(cp3*D);

   //Save the state
		for (s=0;s<S;s++) {P_prev[s]=P[s]; P_m_prev[s]=P_m[s];V_prev[s]=V[s];}
		best_prev=best;

 // Modified bits. Useful for some methods
for (s=0;s<S;s++)
{
  B_mod[s].size=D;
  for (d=0;d<D;d++) B_mod[s].b[d]=0; // Not yet modified
}
   
//-------------------------------------- ITERATIONS
init=0;
loop:
Diam=diameter();
//fprintf(f_trace,"\n%i %i",nb_eval,Diam);
//fprintf(f_trace,"\n%i %i",nb_eval,nb_pivots);
n_iter=n_iter+1;

 if(method==50)
 {



    n_adapt=n_adapt+1;
    if(n_adapt>n_adapt_max)
    {
      n_adapt=0;
      if (error>=error0) k0=k0+1; else {k0=k0-1;if(k0<1) k0=1;}
 //printf("\n %i %f %f",k0, error0,error);     
      error0=error;
    } 
  }


 if (adapt_method==1) 
 {
  if (error>=error_prev)// If there has been no improvement ...
  {

	 n_adapt=n_adapt+1;
	if(n_adapt>n_adapt_max) // ... and if 
	{
    method_rank=(method_rank+1)%strat_number;
	//method_rank=alea_integer(0,strat_number-1); // Test
		method=strategies[method_rank];
		n_adapt=0;
	}
  }
  else n_adapt=0; // TEST
 }


if (adapt_method==2) // Try TWO strategies starting from the same state
					// and keep the best one
{
//fprintf(f_trace,"\n%i %i %f",n_iter,method, error);	
	if(n_iter==n_iter1+adapt_latency*n_adapt_max)
	{
		n_iter0=n_iter;
		error1=error;
   //Save the state
		for (s=0;s<S;s++) {P_prev[s]=P[s]; P_m_prev[s]=P_m[s];V_prev[s]=V[s];}
    //*** (TO DO: add saving for B_mod if "taboo" is use, and restore below)
		best_prev=best;
//fprintf(f_trace,"\n save the state");
	}
	if(n_iter==n_iter0+n_adapt_max)
	{
		error2=error; // Error after iter n_iter-1
//fprintf(f_trace,"\niter %i method %i error2 %f",n_iter,method,error2);

     // Save the final state after having used method rank 0
		for (s=0;s<S;s++) {P_temp[s]=P[s]; P_m_temp[s]=P_m[s];V_temp[s]=V[s];}
		best_temp=best;

		for (s=0;s<S;s++) {P[s]=P_prev[s];P_m[s]=P_m_prev[s];V[s]=V_prev[s];} // Back

		best=best_prev;

		method_rank=(method_rank+1)%2; // Switch method
		method=strategies[method_rank];
		init_links=0; // Keep the same links
	}
	if(n_iter==n_iter0+2*n_adapt_max) 
	{
		error3=error; // Error after back and having tried method rank 0
//fprintf(f_trace,"\niter %i method %i error3 %f",n_iter,method,error3);
		if(error3>error2) // If no improvement
		{			
			method_rank=(method_rank+1)%2; // Switch (i.e. back to the previous method)
			method=strategies[method_rank];
			init_links=0; // Keep the same links

			for (s=0;s<S;s++) {P[s]=P_temp[s];P_m[s]=P_m_temp[s];V[s]=V_temp[s];} // Back
			best=best_temp;
		}

//printf("\n error2 %f error3 %f n_iter %i method %i",error2,error3,n_iter,method);
		n_iter1=n_iter;
		adapt_latency=adapt_latency+1;
	}
 // if(method==100) K=35; else K=alea_integer(2,3);  //TEST
}

//printf("\n method %i for iter %i, error %f",method,n_iter, error);

if(adapt_K==1)

{
	K_choice[K-1][1]=K_choice[K-1][1]+1; // Number of times K has been used

	//cp1=(error_prev2-error_prev)/error_prev2;



	cp2=(error_prev-error) /error_prev;
	//cp2=cp2*cp2;
	cp2=log(1+cp2);
	K_choice[K-1][0]=K_choice[K-1][0]+cp2; // Cumulated improvement

//if (cp2>0) printf("\n %f", cp2);
	n_adaptK=n_adaptK+1; // Let some time to the current K

	if(n_adaptK>n_adapt_max)
	{
		if (K_choice[0][1]*K_choice[2][1]==0)
		{
			K=4-K;// Switch. 
		}
		else
		{
			cp0=K_choice[0][0]/K_choice[0][1];
			cp1=K_choice[2][0]/K_choice[2][1];
			pr=alea(0, cp0+cp1);
			if (pr<cp0) K=1; else K=3;
//printf("\n %.0f %.0f %f %f pr %f =>K %i", K_choice[0][1],K_choice[2][1],cp0,cp1,pr,K);
//printf("\n %i",K);
		}
		n_adaptK=0;
	}
}



if(init_links==1 && !fixed_link) // Who informs who, at random
{
  for (s=0;s<S;s++) //
	{
		for (m=0;m<S;m++) LINKS[m][s]=0;
		LINKS[s][s]=1; // However, each particle informs itself
	}

  // Other links
	  for (m=0;m<S;m++)
	  {
		s1=0;s2=S-1;

		if (K>1)
			for (k=0;k<K;k++)
			{
			s=alea_integer(s1,s2);
			LINKS[m][s]=1;
		//	LINKS[s][m]=1;
			}
	  }
}
if (fixed_link && nb_eval<2*S) // Constant circular neighbourhood
{
	k1=(K-1)/2;
	if(K%2==0) k2=K/2; else k2=k1;
	for (s=0;s<S;s++) //
	{
		for (m=0;m<S;m++) LINKS[m][s]=0;
		for (m=s;m<=s+k2;m++) 
		{i=m%S; LINKS[i][s]=1;}
		for (m=s;m>=s-k1;m=m-1) 
		{i=(S+m)%S; LINKS[i][s]=1;}
	}
}

if(print_level>1) // Display links

{
   printf("\nLINKS");
	for (s=0;s<S;s++) //
	{ printf("\n");
		for (m=0;m<S;m=m++) printf("%i ",LINKS[m][s]);
	}
}


// The swarm MOVES

for (s=0;s<S;s++) pivots[s]=0;

for (s=0;s<S;s++)  // Best informants
{
	g=best_info(s,option_best);
  pivots[g]=1;
    // Note. TO DO. Memorise g in a table and use it below
 }
 // Number of pivots . The mean number is S/(K-1)
 nb_pivots=0; for (s=0;s<S;s++) nb_pivots=nb_pivots+ pivots[s];
 // fprintf(f_trace,"\n %i",nb_pivots);
  
for (s=0;s<S;s++)  // For each particle ...
{  
	g=best_info(s,option_best);
 
switch (method)
{
	
	case 0: // Magic PSO
 /*
	if (S==1) 
		cp2=alea(-1,1);// cp2=alea(0,1); // Test for S=1
	else
	{  
		if(deter==0) // Random choices in pure binary algo
		{  
			i=alea_integer(0,2);
			//i=alea_integer(0,3);
		if (i==0){cp2=1;cp3=1;}   // Towards x_p and x_g
		if (i==1) {cp2=1;cp3=-1;} // Towards x_p away from x_g
		if (i==2) {cp2=-1;cp3=1;} // Away from x_p and towards x_g
		if (i==3) {cp2=-1;cp3=-1;} // Away from x_p and away from x_g
		}

	} 
 */
  cp2=2*alea_integer(0,1)-1;   cp3=2*alea_integer(0,1)-1;
	for (d=0;d<D;d++)
	{
		bit1=P[s].x.b[d]; bit2=P_m[s].x.b[d];

		if (S==1) 
			V[s].b[d]=(int)(V[s].b[d] + 2*cp2*(bit2-bit1)); // Test for S=1
		else
		{
			bit3=P_m[g].x.b[d];

			if(deter==-1) // more randomness in pure binary algo

			{
			i=alea_integer(0,2);
			// i=alea_integer(0,3);
			if (i==0){cp2=1;cp3=1;}   // Towards x_p and x_g
			if (i==1) {cp2=1;cp3=-1;} // Towards x_p away from x_g
			if (i==2) {cp2=-1;cp3=1;} // Away from x_p and towards x_g
			if (i==3) {cp2=-1;cp3=-1;} // Away from x_p and away from x_g
			}
 /*
   if ( V[s].b[d]<-1) // Useful only when using several strategies
   {
    V[s].b[d]=-1; 
    }
    else if (V[s].b[d]>1) V[s].b[d]=1; else V[s].b[d]=0;
 */   
		V[s].b[d]=(int)(V[s].b[d] + cp2*(bit2-bit1) +cp3*(bit3-bit1)); // in [-3...3]
		P[s].x.b[d]=P[s].x.b[d]+V[s].b[d];  // in {-3...4}     
		P[s].x.b[d]=(4+(int)P[s].x.b[d])%2;  // new position in {0,1}
//		V[s].b[d]= V[s].b[d]-bit1+  P[s].x.b[d];
 //printf("\n %f",   V[s].b[d]);
		V[s].b[d]=(int)(3+V[s].b[d])%3-1; // new velocity in {-1,0,1}
 /*       // TEST
    if (V[s].b[d]<-1 || V[s].b[d] >1)
    {printf("\nERROR 0: %f",  V[s].b[d]);
      printf("\n %i %i %i",bit1,bit2,bit3);
      printf("\n %f %f",cp2,cp3);
  scanf("%i",&i);
    }
 */
		}
	}
	break;

case 1:

   // Method 1 -----------------------------------
cp0=0.5;  // Random move, relative probability
if (deter!=-1)
{  
  cp1=cp0;

  cp2=(double)nb_eval*K/(double)S;
  cp2=1-cp1/log(cp2);
  cp3=cp2;
  cp=cp1+cp2+cp3;
}

	for (d=0;d<D;d++)

	{
    if (deter==-1)
    {

      cp1=alea(cp0/2,cp0);
      cp2=(double)nb_eval*K/(double)S;
      cp2=1-cp1/log(cp2);
      cp3=cp2;
      cp=cp1+cp2+cp3;
    }
    bit1=alea_integer(0,1);
    bit2=P_m[s].x.b[d];   bit3=P_m[g].x.b[d];
    
    pr=alea(0,cp);
    if (pr<=cp1) bit= bit1;

    else 
      if (pr<=cp1+cp2) bit=bit2; else bit=bit3;


     P[s].x.b[d]=bit; 
    //V[s].b[d]=  P[s].x.b[d]-bit0; // -1 or 0 or 1. 
	}
break;


case 2:
// Method 2 ---------------------------
cp0=2/3.;

cp0=0.8;
 cp1=cp0;

  cp2=(double)nb_eval*K/(double)S;
  cp2=1-cp1/log(cp2);

  cp3=cp2;

	for (d=0;d<D;d++)
	{
    if (deter==-1)
    {
      cp1=alea(cp0/2,cp0);
      cp2=(double)nb_eval*K/(double)S;
    cp2=1-cp1/log(cp2);
    cp3=cp2;
    }
    
    bit1=alea_integer(0,1);  bit2=P_m[s].x.b[d];   bit3=P_m[g].x.b[d];
    pr=alea(0,1); if (pr<cp1) choice[0]= bit1; else choice[0]=1-bit1;
    pr=alea(0,1); if (pr<cp2) choice[1]=bit2; else choice[1]=1-bit2;
    pr=alea(0,1); if (pr<cp3) choice[2]=bit3; else choice[2]=1-bit3;

    //Majority
    bit=choice[0]+choice[1]+choice[2];
    if(bit<2)  P[s].x.b[d]=0; else P[s].x.b[d]=1;

	}
break;

case 3:
// Method 3 -----------------------
	cp1=0.5;
	cp=(double)nb_eval*K/(double)S;
	cp=cp1/log(cp);


	for (d=0;d<D;d++)
	{
		bit0=0; bit1=0;

		for (s1=0;s1<S;s1++) // Ask all informers
		{
			if(LINKS[s1][s]==0) continue;
				bit=P_m[s1].x.b[d];

				if (bit==0) bit0=bit0+1; else bit1=bit1+1;		
		}

		pr=alea(0,bit0+bit1);
		if(pr<=bit0) bit=0; else bit=1; // Probabilistic majority...

		pr=alea(0,1); if(pr<=cp) bit=1-bit;// ... not sure 
		P[s].x.b[d]=bit;
	}

break;

case 4:
// Method 4 ------------------------------------------
  cp1=0.5;
  cp2=(double)nb_eval*K/(double)S;
  cp2=1-cp1/log(cp2);
  cp=cp1+K*cp2;

	for (d=0;d<D;d++)
	{
		bit0=0; bit1=0;
		for (s1=0;s1<S;s1++) // Ask all informers
		{
			if(LINKS[s1][s]==0) continue;

				bit=P_m[s1].x.b[d];
				if (bit==0) bit0=bit0+1; else bit1=bit1+1;
		}
     // Weighted majority
		if(bit0>bit1) bit=0;else if (bit1>bit0) bit=1; else bit=alea_integer(0,1);

    pr=alea(0,cp); if(pr<cp1) bit=alea_integer(0,1);
		P[s].x.b[d]=bit;
	}
break;


case 5:
// Method 5---------------------------------

  //cp1: Confidence coefficient for velocity
  //cp2: Confidence coefficient for p
  //cp3: Confidence coefficient for g

cp0=1-0.5*log((double)D)/log((double)S);
 if (K>1) cp0=cp0*log((double)K);


 z=(double)nb_eval*K/(double)S;
 if (z>1) z=log(z);


 cp1=cp0/z; // Decreasing. Equal to initial cp0 just after initialisation

 cp2=cp1; // Decreasing
 if(K>2)
  cp3=cp2/log(K-1);// Decreasing. Slightly smaller than cp2
 else cp3=cp2;

 

   n1=(int)(cp1*D);   // Decreasing.
   n2=(int)(cp2*D);// Decreasing. The more you trust p, the less you update it
   n3=(int)(cp3*D); //Decreasing.  The more you trust g, the less you update it

  Pp=P_m[s];Gp=P_m[g];

  for (d=0;d<D;d++)   Vp.x.b[d]=-1;  // Arbitrary value: means "non assigned"

  for (k=0;k<n1;k++) // Velocity. At most n1 non assigned components
  {
    d=alea_integer(0,D-1); Vp.x.b[d]=alea_integer(0,1);
  }
  for (k=0;k<n2;k++)  // Around p. At most n2 modified components
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
    for (k=0;k<n3;k++) // Around g. At most n3 modified components
  {
    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  } 
 
	for (d=0;d<D;d++) // Choose each component
	{
		bit1=Vp.x.b[d];  bit2=Pp.x.b[d];   bit3=Gp.x.b[d];

		//Majority choice
		if(bit1>=0) // If assigned velocity
		{ bit=bit1+bit2+bit3;  
		if(bit<2)  P[s].x.b[d]=0; else P[s].x.b[d]=1;
		}
		else // If non assigned velocity
		{
			 bit=bit2+bit3;
			if(bit==2) P[s].x.b[d]=1; else
			{if(bit==0) P[s].x.b[d]=0; else P[s].x.b[d]=alea_integer(0,1);}
		}
	}
break;


case 6:
// Method 6---------------------------------
  //cp2: Confidence coefficient for p
  //cp3: Confidence coefficient for g

cp0=1-0.5*log((double)D)/log((double)S);
 cp0=cp0+((double)K/(double)S)*log((double)D)/log((double)S);
if (K>1) cp0=cp0*log((double)K);

 z=(double)nb_eval*K/(double)S;
 if (z>1) z=log(z);

 cp1=cp0/z;
 cp2=cp1;


 if(K>2)
  cp3=cp1/log(K-1);// Decreasing. 
 else cp3=cp1;

   n2=(int)cp2*D;// Decreasing. The more you trust p, the less you update it
   n3=(int)cp3*D; //Decreasing.  The more you trust g, the less you update it

  Pp=P_m[s];Gp=P_m[g];

  for (k=0;k<n2;k++)  // Around p. At most n2 modified components
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
    for (k=0;k<n3;k++) // Around g. At most n3 modified components
  {
    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  } 
  
	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=Pp.x.b[d];   bit3=Gp.x.b[d];


		//Majority choice
	//	if(bit2==bit3 && alea(0,1)<cp1) P[s].x.b[d]=bit2;
    if(bit2==bit3) P[s].x.b[d]=bit2;  
		else P[s].x.b[d]=alea_integer(0,1);
	}
break;

case 7:
// Method 7---------------------------------
// Like 5, but coefficients are modified only if there has been


// some improvement

  //cp2: Confidence coefficient for p
  //cp3: Confidence coefficient for g

if (error<error_prev) goto skip7;

	z7=z7+(double)K/(double)S;
	z=log(z7);
	cp1=cp0/z;

	cp2=cp1;

	if(K>2)

	cp3=cp1/log(K-1);// Decreasing. 
	else cp3=cp1;

   n2=(int)cp2*D;// Decreasing. The more you trust p, the less you update it
   n3=(int)cp3*D; //Decreasing.  The more you trust g, the less you update it
skip7:

  Pp=P_m[s];Gp=P_m[g];

  for (k=0;k<n2;k++)  // Around p. At most n2 modified components
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
    for (k=0;k<n3;k++) // Around g. At most n3 modified components
  {
    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  } 
  
	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=Pp.x.b[d];   bit3=Gp.x.b[d];

		//Majority choice
		if(bit2==bit3) P[s].x.b[d]=bit2;	
		else P[s].x.b[d]=alea_integer(0,1);
	}
break;

case 71:
// Like 7, but with probabilistic majority choice
if (error<error_prev) goto skip71;

	z7=z7+(double)K/(double)S;
	z=log(z7);
	cp1=cp0/z;
	cp2=cp1;

	if(K>2)

	cp3=cp1/log(K-1);// Decreasing. 
	else cp3=cp1;

   n2=(int)cp2*D;// Decreasing. The more you trust p, the less you update it
   n3=(int)cp3*D; //Decreasing.  The more you trust g, the less you update it
skip71:

  Pp=P_m[s];Gp=P_m[g];

  for (k=0;k<n2;k++)  // Around p. At most n2 modified components
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
    for (k=0;k<n3;k++) // Around g. At most n3 modified components
  {

    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  }
  
	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=Pp.x.b[d];   bit3=Gp.x.b[d];
		bit=bit2+bit3;
		//Probabilistic majority choice
		P[s].x.b[d]=alea_integer(0,1);
		pr=alea(0,1);
		// You may have 00 just by chance.
		if (bit==0 && pr>0.025*cp2*cp3) {P[s].x.b[d]=0;continue;}
		// You may have 11 just by chance
		if (bit==2 && pr>0.025*cp2*cp3) {P[s].x.b[d]=1;continue;}

//continue;
		// method 71a: to be more rigourous, you may also take into account
		// 01 and 01: 

		if (bit2==0 && bit3==1 && pr<0.5*cp2*(1-cp3)) {P[s].x.b[d]=1;continue;}
		if (bit2==0 && bit3==1 && pr<0.5*(1-cp2)*cp3) {P[s].x.b[d]=0;continue;}
		if (bit2==1 && bit3==0 && pr<0.5*(1-cp2)*cp3) {P[s].x.b[d]=1;continue;}
		if (bit2==1 && bit3==0 && pr<0.5*cp2*(1-cp3)) {P[s].x.b[d]=0;continue;}
		
	}
break;

case 72: // Mixing 7 and 11. Constant coefficients

    dist=(int)log(D);
	  
	Gp=P_m[g];
	r=alea(1,dist);

    for (k=0;k<r;k++) // Around g.
  {

    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d]; 
  } 
  
	if(s==g) { P[s]=Gp; break;}
 
	Pp=P_m[s];	
	r=alea(1,dist);
  for (k=0;k<r;k++)  // Around p.
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }

	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=Pp.x.b[d];   bit3=Gp.x.b[d];

		//Majority choice
		if(bit2==bit3) P[s].x.b[d]=bit2;	
		else 
		{	
			P[s].x.b[d]=alea_integer(0,1);
	/*		
		// "a" variant. Trust g a bit more
			pr=alea(0,1+1/log(K-1)); 
			if (pr<=1) P[s].x.b[d]=bit2; else P[s].x.b[d]=bit3;
	*/	
		}
	}
break;

case 73:
if (error<error_prev) goto skip73;

	z7=z7+(double)K/(double)S;
	z=log(z7);
	cp1=cp0/z;
	cp2=cp1;

	if(K>2)

	cp3=cp1/log(K-1);// Decreasing. 
	else cp3=cp1;

skip73:

  Pp=P_m[s];Gp=P_m[g];
	dist=(int)cp2*D;
	r=alea(1,dist);
  for (k=0;k<r;k++)  // Around p.
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
	dist=(int)cp3*D;
	r=alea(1,dist);
    for (k=0;k<r;k++) // Around g. 
  {
    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  } 
  
	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=Pp.x.b[d];   bit3=Gp.x.b[d];

		//Majority choice
			if(bit2==bit3) P[s].x.b[d]=bit2;	
		else   // P[s].x.b[d]=alea_integer(0,1);
		{	
			pr=alea(0,1+log(K-1));
			if (pr<=1) P[s].x.b[d]=bit2; else P[s].x.b[d]=bit3;
		}
		
	}
break;

case 8:
// Method 8---------------------------------

  Pp=P_m[s];Gp=P_m[g];


	for (d=0;d<D;d++) // Choose each component
	{
     bit2=Pp.x.b[d];   bit3=Gp.x.b[d];

    //Majority choice
		if(bit2==bit3) P[s].x.b[d]=bit2; else P[s].x.b[d]=alea_integer(0,1);
	}
break;

case 9:
// Method 9 ------------------------------------------------

//pr=1.0/pow(2,D);

c2=(double)nb_eval/(double)S;
c3=K*c2;
//cp2=pow(1-pr,c2);
//cp3=pow(1-pr,c3);
cp2=c2*log(1-pr);  // log(Failure proba)
cp3=c3*log(1-pr);

//printf("\n %f %f ",cp2,cp3);

if(s!=g)
{
	for (d=0;d<D;d++) // Choose each component
	{
		 bit2=P_m[s].x.b[d];   bit3=P_m[g].x.b[d];

     if (log(alea(0,1))>cp2) bit2=1-bit2;

     if(log(alea(0,1))>cp3) bit3=1-bit3;
		//Majority choice
    if(bit2==bit3) P[s].x.b[d]=bit2;
		else P[s].x.b[d]=alea_integer(0,1);
	}
}
else
{	
	for (d=0;d<D;d++) // Choose each component
	{
		bit3=P_m[g].x.b[d];
		if (alea(0,1)>cp3) bit3=1-bit3;
		P[s].x.b[d]=bit3;
	}
}
break;

case 10:
// Method 10 ------------------------------------------ 
cp0=0.5;
 cp1=cp0;
  cp2=(double)nb_eval*K/(double)S;
  cp2=1-cp1/log(cp2);
  cp3=cp2;

	for (d=0;d<D;d++)
	{
    if (deter==-1)
    {
      cp1=alea(cp0/2,cp0);
      cp2=(double)nb_eval*K/(double)S;
    cp2=1-cp1/log(cp2);
    cp3=cp2;
    }

    bit1=P[s].x.b[d];  bit2=P_m[s].x.b[d];   bit3=P_m[g].x.b[d];
    pr=alea(0,1); if (pr<cp1) choice[0]= bit1; else choice[0]=1-bit1;
    pr=alea(0,1); if (pr<cp2) choice[1]=bit2; else choice[1]=1-bit2;
    pr=alea(0,1); if (pr<cp3) choice[2]=bit3; else choice[2]=1-bit3;

    //Majority
    bit=choice[0]+choice[1]+choice[2];
    if(bit<2)  P[s].x.b[d]=0; else P[s].x.b[d]=1;
    V[s].b[d]=  P[s].x.b[d]-bit1; // -1 or 0 or 1. For information
	}

  break;

case 11:  // Pivot method, from an idea of P. Serra, 1997
          // (seriously) modified and adapted to binary problems by M. Clerc, 2004
case11:
	// Works pretty well on some problems .. and pretty bad on some others
  // Note that there is no velocity, and the current position is not used
  // It works better when information links are modified at random
  // if there has been no improvement after the previous iteration
  
 P[s]=P_m[g]; // Best previous position g of the best informant

 // (for D>=2) Simplified formula.  Note that "dist" is an integer
 //dist=log(D);  if (dist<3) dist=3;  // 11
 //dist=D*log(D)/S;// 11e
 //dist=1+(K-1)*Diam/(double)S; //if (dist<2)dist=2;   // 11f
 //printf("\n %i %i",nb_pivots, Diam);
 dist=(int)(1+Diam/(double)nb_pivots); if (dist<3) dist=3; //* 11g

//if(dist>3)printf("\n %i",dist); 
 if (dist>D) dist=D;  // For D=1.

  r=alea(1,dist); //11. Note that r is never equal to dist
 // r=1+alea_integer(1,dist-1);  //11b. Theoretically equivalent
 // r=1+alea_integer(1,dist); //11c. A bit larger
 // r=alea_integer(1,dist-1);  // 11d.   More uniform , on dist k values
                              // Sometimes far better but also sometimes far worse
  
 for (d=0;d<D;d++) distrib[d]=0; // Just for some stats
 
 for (k=0;k<r;k++) // Define a position around g
 //    Note that even if there at least two k values
 //  there is a small probability (1/D^k)
 // that just one bit is modified
  {
    d=alea_integer(0,D-1);
    P[s].x.b[d]=1- P[s].x.b[d];
  distrib[d]=1;
  }
   i=0;    for (d=0;d<D;d++) i=i+distrib[d];
   k_tot[i]=k_tot[i]+1; // Just for some statistics
 
 break;

case 12:  // Pivot method, variant "independent particles"
		// K is not taken into account
 Gp=P_m[s];

 dist=(int)log(D); 

  r=alea(1,dist);
 //r=1+alea_normal(0,dist);// Test

 for (k=0;k<r;k++)
  {
    d=alea_integer(0,D-1);
    Gp.x.b[d]=1- Gp.x.b[d]; // Around s
  }
    P[s]=Gp;
 break;


 case 13:  // Pivot method, variant "one random informer"
		// K is not taken into account
d=alea_integer(0,D-1);
if (total_error(P_m[d].f)<total_error(P_m[s].f)) g=d; else g=s;
 Gp=P_m[g];

  dist=(int)log(D); // c

  r=alea(1,dist);
 //r=1+alea_normal(0,dist);// Test

 for (k=0;k<r;k++)
  {
    d=alea_integer(0,D-1);
    Gp.x.b[d]=1- Gp.x.b[d]; // Around s
  }
    P[s]=Gp;
 break;


case 14: 

case14:
/*
 for (i=0;i<S;i++)  // For each particle ...
{  
    g_[i]=best_info(i,option_best); // .. find the best informer
} 
	
// Find the nearest one
 z=D+1;
 for (i=0;i<S;i++)
 {
	m=g_[i];
	if(m==g) continue;
		r=distance(g,m,-2);
		if(r<=0) continue;
			if(r<z) z=r;
 }
*/
/*

 z=D+1; // Find the nearest informer
 for (m=0;m<S;m++)
 {
	if(m==g) continue;
		r=distance(g,m,-2);
		if(r<=0) continue;
			if(r<z) z=r;
 }

*/

//----------- 
	// Good for M�hlenbein, Multimodal
	// Bad for Zebra3
 Gp=P_m[g];

  // A "center" between s and g
 for (d=0;d<D;d++) 
 {
	 if (P_m[s].x.b[d]!=P_m[g].x.b[d] && alea(0,1)<0.5) Gp.x.b[d]=P_m[s].x.b[d];
 }
z=distance(s,g,-1)+1;
r=alea(1,z);

 for (k=0;k<r;k++)
  {

    d=alea_integer(0,D-1);
    Gp.x.b[d]=1- Gp.x.b[d]; // Around g
  }
    P[s]=Gp;
 break;

case 15:
 if (error<error_prev) goto case14; else goto case11;


case 16: // Like 11, but with "taboo"
Gp=P_m[g];

  dist=(int)log(D);  // We suppose here D>=2
  r=alea(1,dist);

 for (k=0;k<r;k++)
  {
    d=alea_integer(0,D-1);
    if(B_mod[s].b[d]<1) // "Taboo". Modify only if it has not just been modified
      {
        Gp.x.b[d]=1- Gp.x.b[d]; // Around g
      }
  }
 
    for (d=0;d<D;d++)

    {
      if (Gp.x.b[d]!=P[s].x.b[d])  B_mod[s].b[d]=1; // Has been modified
          else B_mod[s].b[d]=0; // Has not been modified
       P[s].x.b[d] =Gp.x.b[d];
      }

 
 break;

 case 17: // Like 16, but with condition "the particle has just improved its position"
Gp=P_m[g];
  dist=(int)log(D);  // We suppose here D>=2
  r=alea(1,dist);
  error1=total_error(P[s].f);
  error2=total_error(P_m[s].f);
   i=error1<=error2; // Improvement test
   
 for (k=0;k<r;k++)
  {
    d=alea_integer(0,D-1);

    if((i && B_mod[s].b[d]<1) || !i)
      {
        Gp.x.b[d]=1- Gp.x.b[d]; // Around g
      }
  }

    for (d=0;d<D;d++)
    {
      if (Gp.x.b[d]!=P[s].x.b[d])  B_mod[s].b[d]=1; // Has been modified
          else    B_mod[s].b[d]=0; // Has not been modified
       P[s].x.b[d] =Gp.x.b[d];
      }
 break;



case 19: // Like 11, but with "taboo"
  //dist=log(D);  // We suppose here D>=2
   dist=(int)(1+Diam/(double)nb_pivots); if (dist<3) dist=3;
  r=alea(1,dist);
 P[s]=P_m[g];
  error1=total_error(P[g].f);
  error2=total_error(P_m[g].f);

  i=error1<=error2; // g has just improved its position
  
 for (k=0;k<r;k++)   // Around g
  {
    d=alea_integer(0,D-1);
    if(!i || (i && B_mod[s].b[d]<1)) // "Taboo"
      {
        P[s].x.b[d]=1- P[s].x.b[d]; 
        B_mod[s].b[d]=1;

      }
    else
    {
      B_mod[s].b[d]=0;
    }

  }

 break;
 
case 99:   // TEST  . 7 with taboo

if (error<error_prev) goto skip99;


	z7=z7+(double)K/(double)S;

	z=log(z7);
	cp1=cp0/z;
	cp2=cp1;


	if(K>2)

	cp3=cp1/log(K-1);// Decreasing.
	else cp3=cp1;

   n2=(int)cp2*D;// Decreasing. The more you trust p, the less you update it
   n3=(int)cp3*D; //Decreasing.  The more you trust g, the less you update it
skip99:

  Pp=P_m[s];Gp=P_m[g];

  for (k=0;k<n2;k++)  // Around p. At most n2 modified components
  {
    d=alea_integer(0,D-1); Pp.x.b[d]=1- Pp.x.b[d];
  }
    for (k=0;k<n3;k++) // Around g. At most n3 modified components
  {
    d=alea_integer(0,D-1); Gp.x.b[d]=1- Gp.x.b[d];
  }

	for (d=0;d<D;d++) // Choose each component
	{
     if(B_mod[s].b[d]<1) // Taboo

     {
         bit1=P[s].x.b[d];
        bit2=Pp.x.b[d];   bit3=Gp.x.b[d];


		    //Majority choice
		    if(bit2==bit3) P[s].x.b[d]=bit2;
		    else P[s].x.b[d]=alea_integer(0,1);

        if (P[s].x.b[d]!=bit1) B_mod[s].b[d]=1;
      }
      else B_mod[s].b[d]=0;
	}
break;
  

//=========== Methods a priori designed for constant circular neighbourhood
// But, of course you may try them with random link (cf fixed_link parameter)
case 100: //  KE-C Method (Jim Kennedy and Russ Eberhart method, improved by Maurice Clerc)
	cp0=1;
	cp1=2.0;
	cp2=cp1;

	for (d=0;d<D;d++)
    {
        cp3=cp0*V[s].b[d] + alea(0,cp1)*(P_m[s].x.b[d]-P[s].x.b[d]);  // cp0*v + rand(0,cp1)*(p-x)
        cp3=cp3+alea(0,cp2)*(P_m[g].x.b[d]-P[s].x.b[d]); //  + rand(0,cp2)*(g-x)
        //Note that g is either the local best or the global best
        // depending on the neighbourhood size
        
        V[s].b[d]=(int)cp3; // New integer velocity
        cp3=  P[s].x.b[d]+ V[s].b[d]; // Add velocity
        cp3=1/(1+exp(-cp3)); // Keep position in ]0,1[

        pr=alea(0,1); // Probabilistic binary choice
        if (pr<cp3) P[s].x.b[d]=1; else P[s].x.b[d]=0;
    }
    break;

 case 101: // TL Method (Tasgetiren and Liang. 2004)
          // Suggested S=2*D, K=S
  cp0=4.0; // Max velocity (absolute value)
	cp1=2.0;
	cp2=cp1;


	for (d=0;d<D;d++)
    {
        cp3=V[s].b[d]+ alea(0,cp1)*(P_m[s].x.b[d]-P[s].x.b[d]);  // cp0*v + rand(0,cp1)*(p-x)
        cp3=cp3+alea(0,cp2)*(P_m[g].x.b[d]-P[s].x.b[d]); //  + rand(0,cp2)*(g-x)
                                      //Note that g is either the local best or the global best
                                      // depending on the neighbourhood size
         if(cp3>cp0) cp3=cp0; else if(cp3<-cp0) cp3=-cp0;  // Take vmax into account
          V[s].b[d]=(int)cp3 ;   // New velocity
          
         // Probabilistic binary choice for position 
        cp3=1/(1+exp(-cp3)); // Keep in ]0,1[
        pr=alea(0,1); 
        if (pr<cp3) P[s].x.b[d]=1; else P[s].x.b[d]=0;
    }
    break;
    
case 111:  // Pivot method -----------------------------------------

	// Works pretty well on some problems .. and pretty bad on some others
 Gp=P_m[g]; // Best previous position g of the best informer

  dist=(int)log(D);  // We suppose here D>=2

  r=alea(1,dist);
 //r=1+alea_normal(0,dist);// Test

 for (k=0;k<r;k++)
  {
    d=alea_integer(0,D-1);
    Gp.x.b[d]=1- Gp.x.b[d]; // Around g
  }
    P[s]=Gp;
 break;

 default:
 printf("\n SORRY, no method %i",method);
 break;
} // End of switch (method)


	// ... evaluate the new position
	for (i=0;i<n_f;i++) P[s].f.fi[i]=fabs(perf(s,function[i])-fMin.fi[i]);


switch (tribe)
{
  default: // ... update the best position
	if (best_(P[s].f,P_m[s].f)==1) P_m[s]=P[s];
  break;
  
  case 1: // ... or update the best position of the whole tribe
  for(s1=0;s1<S;s1++)
  {
    	if (best_(P[s].f,P_m[s1].f)==1) P_m[s1]=P[s];
  }
  break;


  case 2:  // update the worst
  i=s; // rank of the worst
  for(s1=0;s1<S;s1++)
  {
     if(s1==s) continue;
    if(LINKS[s1][s]==0) continue;
    if (best_(P_m[s1].f,P_m[i].f)==1) continue;
      i=s1; 
  }
// if(i!=s) printf("\n %i modify %i",s,i);
    if (best_(P[s].f,P_m[i].f)==1) P_m[i]=P[s];
  break;
 } // end switch (tribe)
 
	// ... update the best of the bests
	if (best_(P_m[s].f,best.f)) best=P_m[s];

// End of evaluation

  if(print_level>0)
{

  printf("\ncurrent ");
   for (d=0;d<D;d++) printf("%1i",(int)P[s].x.b[d]);
  for (i=0;i<n_f;i++) printf(" %f",P[s].f.fi[i]);
  printf("\n  best ");
     for (d=0;d<D;d++) printf("%1i",(int)P_m[s].x.b[d]);
  for (i=0;i<n_f;i++) printf(" %f",P_m[s].f.fi[i]);
  printf("\n    best best ");
  for (d=0;d<D;d++) printf("%1i",(int)best.x.b[d]);
  for (i=0;i<n_f;i++) printf(" %f",best.f.fi[i]);
}


}


// Check if finished
error_prev2=error_prev;
error_prev=error;
error=total_error(best.f);
//printf("\n error %f",error);

/*
// If no improvement, information links will be reinitialized
if(error<error_prev)
{
  init_links=0;
}
else

{
  init_links=1;
}
*/
	
init_links=1;

//fprintf(f_trace,"\n%i %f",nb_eval,save_time);

if (error>eps && nb_eval<eval_max) goto loop;
if (error>eps) n_failure=n_failure+1;

// Result display-save
printf("\n\nExec %i Eval= %i ",n_exec,nb_eval);
printf("\n");
for (d=0;d<D;d++) printf("%1i",(int)best.x.b[d]);
for (i=0;i<n_f;i++) printf(" %.4f",best.f.fi[i]);

if (n_exec==1) min=error; else if(error<min) min=error;

if (function[0]==11) // Special save for this case
{
	DD=D/Sp;
	for (i=0;i<Sp;i++)
	{
		printf("\n s%i ",i+1);
		for (d=0;d<DD;d++) printf("%i",best.x.b[d+DD*i]);

		for (d=0;d<DD;d++) fprintf(f_trace,"%i ",best.x.b[d+DD*i]);
		fprintf(f_trace,"\n");
	}

	for (i=0;i<Sp-1;i++) // "particle" i
	{
		printf("\n");
		for (j=i+1;j<Sp;j++)  // "particle j
		{
			// Compute the distance
			dist=0;
			for (d=0;d<DD;d++)
			{
				if (best.x.b[d+DD*i]!=best.x.b[d+DD*j]) dist=dist+1;
			}
			printf(" %i",dist);
		}
	}
	fprintf(f_trace,"\n");
}

// Save error for stats
error_[n_exec]=error;
printf("\n %i %f",n_exec,error_[n_exec]);

// Multiobjective result save
if (n_f>1)
{
	for (i=0;i<n_f;i++)
		fprintf(f_multi,"%f ",best.f.fi[i]);
	for (d=0;d<D;d++) fprintf(f_multi," %.4i",(int)best.x.b[d]);
	fprintf(f_multi,"\n");
}


if (adapt_K==1) printf("\n %f %f",K_choice[0][0],K_choice[2][0]);

// Compute and display some statistical information


eval_mean=eval_mean+nb_eval;
eps_mean=eps_mean+error;

if (n_exec<n_exec_max) goto init;

t2=clock();
printf("\n\n Time (clocks) %.0f",t2-t1);
save_time= save_time/eval_mean;
printf("\n Saved time rate %f",save_time);

eval_mean=eval_mean/(double)n_exec;
eps_mean=eps_mean/(double)n_exec;
printf("\n Eval. (mean)= %f",eval_mean);
printf("\n Error (mean) = %f",eps_mean);
success_rate=1- n_failure/(double)n_exec;
printf("\n Success rate = %.6f%%",100*success_rate);
if (n_exec>1) printf("\n Best min value = %f",min);

r=0;
for(k=1;k<=n_exec;k++) 
{r=r+error_[k];
	//printf("\n k %i r %f",k,r);
}

r=r/n_exec; // Mean
pr=0 ;
for(k=1;k<=n_exec;k++) pr=pr+(r-error_[k])*(r-error_[k]);
pr=sqrt(pr/n_exec);
printf("\n Mean %f  StDev %f",r,pr);

fprintf(f_synth,"%f ",eval_mean);
fprintf(f_synth,"%f ",eps_mean);

fprintf(f_synth,"%.6f%% ",100*success_rate);

fprintf(f_synth,"%f ",min);
fprintf(f_synth,"%f ",save_time);
fprintf(f_synth,"\n");


if (eval_max==eval_max1) quality=success_rate*(alpha-pow(alpha,eval_max1-1));
else
	if (eval_max<eval_max2)
		quality=quality+success_rate*(pow(alpha,eval_max-delta_eval_max-1)-pow(alpha,eval_max-1));
	else
		quality=quality+success_rate*(pow(alpha,eval_max2-delta_eval_max-1));
//fprintf(f_synth,"\n %i %f ",eval_max,quality);

} // End loop on eval_max

z=quality/alpha;
printf("\n Quality: %f",z);
fprintf(f_synth,"\nQuality %f ",z);

printf("\n k-distribution stat."); for (k=1;k<=dist+1;k++) printf(" %.0f",k_tot[k]);

end:;
}



//===========================================================
double	alea(double a,double b)
{	// rand number (uniform distribution) in [a b]

double 	r;
 r=(double)rand_kiss()/RAND_MAX_KISS;
// r=(double)rand(); r=r/RAND_MAX;
return a+r*(b-a);
}
//===========================================================
int	alea_integer(int a,int b)
{// Integer rand number in [a b]
int	ir;
double 	r;
r=alea(0,1); ir=(int)(a+r*(b+1-a)); if (ir>b) ir=b;
return ir;
}
//===========================================================
 double alea_normal(double mean, double std_dev)
 {
 /*
     Use the polar form of the Box-Muller transformation to obtain
     a pseudo random number from a Gaussian distribution
	 Note that it is always >= mean

 */
        double x1, x2, w, y1;
        //double  y2;

         do {
                 x1 = 2.0 * alea(0,1) - 1.0;

                 x2 = 2.0 * alea(0,1) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0);

         w = sqrt( -2.0 * log( w ) / w );
         y1 = x1 * w;

        // y2 = x2 * w;
          y1=y1*std_dev+mean;
         return y1;
  }



//===========================================================

int	best_(struct f fx,struct f fy)
{
	/*
	Compare two position evaluations. In case of multiobjective
	it is the Pareto dominance
	Return 1 is fx better than fy, 0 else
	  */

	int	i;
	for(i=0;i<fx.size;i++)
 {
    if (fx.fi[i]>fy.fi[i]) return 0;
 }
 // All fx values are <=

for(i=0;i<fx.size;i++)
 {
      if (fx.fi[i]<fy.fi[i]) return 1;  // At least one is <
 }
return 0;    //You may also return 2, if this particular case
            // is interesting
}


//===========================================================
int	best_info(int s,int option)
{
/*
Return the rank of the best informer
option	= 0 => classical (really the best)
		= 1 => using pseudo-gradient
*/
double	delta1,delta2;
double	dist;
double	error_1,error_2;
int	g;

double	grad;
int	m;

struct f	min_f;
int	s2;

switch (option)
{

case 1: // Pseudo gradient
	g=s; // If nothing better, the best informer is
			// the particle itself
	grad=0;

	error_1=total_error(P_m[s].f);


	for (s2=0;s2<S;s2++)

	{


		if (LINKS[s2][s]==0) continue; // No link


		error_2=total_error(P_m[s2].f);
		delta1=error_1-error_2;
		dist=distance(s,s2,0);

		if (dist<=0) continue;

		delta2=delta1/dist;
		if (delta2>grad)
		{
			g=s2; grad=delta2;
		}
	}
	// if (g==s) goto case0; // Possible variant
	break;


case 0: // Classical option: the "real" best
	for (g=0;g<S;g++) {if (LINKS[g][s]==0) continue; goto next;}
next:
	min_f=P_m[g].f;

	for (m=g+1;m<S;m++)
	{
		if (LINKS[m][s]==0) continue;
			if (best_(P_m[m].f,min_f)==1) {g=m;min_f=P_m[m].f;}
	}
break;
}

return g;
}

//===========================================================

struct bits binary(unsigned int a, int size)
{
/* Compute the bit string (truncated at "size")
 corresponding to the integer a
*/
 int i,max,rank,z1,z2;
 struct bits bit;

 for (i=0;i<size;i++) bit.b[i]=0;

 bit.size=size;
 if(a==0) return bit;
 if(a==1) {bit.b[0]=1; return bit; }
 
 z1=a;
 rank=0;
 max=1;
 loop:
 z2=z1/2;
 if(z2==0) {bit.b[rank]=1; goto end;}

if (rank<=size-1)
  {
    if(2*z2<z1) bit.b[rank]=1;
    z1=z2; rank=rank+1; goto loop;
  }

end:
// for (i=0;i<size;i++) printf("%1i",bit.b[size-i-1]);

return bit;
}

//===========================================================
double combin(int Dim,int k)
{
   double c;
   int d;
  c=Dim;
  for (d=1;d<=k-1;d++) c=c*(Dim-d)/(double)(d+1);
  return c;
}
//===========================================================
int diameter()   // Swarm diameter
{
 int diam,dist;
 int s1,s2;
 diam=0;
 for (s1=0;s1<S-1; s1++)
 for (s2=s1+1;s2<S;s2++)
 {
   // dist=distance(s1,s2,-1); // Memory swarm
   dist=(int)distance(s1,s2,-2); // Current swarm

   if(dist>diam) diam=dist;
 }
  return diam;
}
//===========================================================
double	distance(int s1,int s2, int option)
{
// Distance between two positions
// option = 0 => s1 and s2 are ranks of best memorized positions
// option = 1 => s1 is the rank of a best memorized position
//               and s2 the rank of a current position
// option <0 => Hamming distance
// WARNING. Depending on your computer, it may not work
// for high dimension
  int	d;
  double dist;
  int n1;
  
  if(option<0) goto hamming;

  dist=0;
  for (d=0;d<P[s1].x.size;d++)
  {
    if(option==0) n1=abs(P_m[s1].x.b[d]- P_m[s2].x.b[d]);
    else n1=abs(P_m[s1].x.b[d]- P[s2].x.b[d]);
    if (n1>0)
    {
      dist=dist+pow(2.0,d);

    }
  }
  if (dist>mean_val) dist=2*mean_val-dist; // Cyclic distance
  return dist;

hamming:
 switch(option)

 {
case -3:
  dist=0;
  for (d=0;d<P_m[s1].x.size;d++)
  {

    n1=abs(P[s1].x.b[d]- P_m[s2].x.b[d]);
	dist=dist+n1;
  }

  return dist;


 default: // Between  current positions
  dist=0;
  for (d=0;d<P[s1].x.size;d++)
  {
    n1=abs(P[s1].x.b[d]- P[s2].x.b[d]);
	dist=dist+n1;
  }
  return dist;

 case -1:
  dist=0;
  for (d=0;d<P_m[s1].x.size;d++)
  {
    n1=abs(P_m[s1].x.b[d]- P_m[s2].x.b[d]);
	dist=dist+n1;
  }
  return dist;
 }

}
//===========================================================
void  fundamental_hyp(int function)
{  // Statistically checks the fundamental hypothesis
 // "Nearer is better"
 // NOTE: this program is in fact "independent". No optimisation

   int  bit;
   int d;
  double fx,fy;
  double f[3];
  int m,n;
//  double dis[3];
  double dist,distx,disty;
//  int loop;
 double  n_sample;
 int nn;
 int option;
// int r;
 int radius;
 float rate_sample;
 double rho;
 int s0,s1;
 double  true;
 double v;

//option: 0 => use distance to the solution
           // (check global improvement)
   // 1 => check local improvement
 

 printf("\n Global 0, local 1, local 2, local 3 ?: ");
 scanf("%i", &option);
 
 printf("\n Sampling rate? ");
 scanf("%f", &rate_sample);

  v=pow(pow(2,D),2);
  n_sample=v*rate_sample;
 //printf("\nDimension is %i. Radius? ",D);
 //scanf("%i",&radius);
  P[0].x.size=D;  P[1].x.size=D;P[2].x.size=D;
    true=0;
switch(option)
{

    case 0:  // For global search
  switch(function) // Solutions are supposed to be known
  {
    case 1:
    for (d=0;d<D;d++) P[0].x.b[d]=1;
    break;

    case 2:
    // Note that there are in fact  pow(D/12,2) solutions

    for(n=0;n<=D/12;n=n+1)
    {
      bit=alea_integer(0,1);

      m=12*n;
      for (d=m;d<m+6;d++) P[0].x.b[d]=bit;
      for (d=m+6;d<m+12;d++) P[0].x.b[d]=1-bit;
    }
    break;

  case 3:
  for (d=0;d<D;d++) P[0].x.b[d]=0;
  break;

  case 6:
  for (d=0;d<D;d++) P[0].x.b[d]=0;
  for(m=3;m<=D;m=m+6)
  {
    for (d=m;d<m+3;d++) P[0].x.b[d]=1;
  }
  break;
  

  case 7:
  for (d=0;d<D-1;d++) P[0].x.b[d]=0;
   P[0].x.b[D-1]=1;
  break;

  case 9:
  for (d=0;d<D;d++) P[0].x.b[d]=0;
  break;

  case -1:
  // Solution for Q12: 000010101001
  for (d=0;d<D;d++) P[0].x.b[d]=0;
   P[0].x.b[4]=1;
   P[0].x.b[6]=1;
   P[0].x.b[8]=1;
   P[0].x.b[11]=1; 
  break;

  case 8:
  // Solution for f_landscape_12_20: 011111111110
    for (d=0;d<D;d++) P[0].x.b[d]=1;
   P[0].x.b[0]=0;
   P[0].x.b[11]=0;
  break;
  
  case 99:
  for (d=0;d<D;d++) P[0].x.b[d]=0;
  break;
  
  default:

  printf("\n NOT YET for function %i, SORRY",function) ;
  goto end;
 }  // end switch function
 //---------------
   printf("\n Solution\n");
 for(d=0;d<D;d++) printf("%i", P[0].x.b[d]); printf("\n");
       // Compute "solution" fitness
      f[0]=perf(0,function);

 for (radius=2;radius<=D;radius++)
 {
    nn=0;
 for (n=0;n<n_sample;n++)
 { 

      //Define at random point x
      for (d=0;d<D;d++) P[1].x.b[d]=alea_integer(0,1);
 
      // Define at random y at distance <=radius

 // for (d=0;d<D;d++) P[2].x.b[d]=alea_integer(0,1);
 
      P[2]=P[1];
      dist=alea_integer(1,radius); 
      for (m=0;m<dist;m++)
      {
          d=alea_integer(0,D-1);
          P[2].x.b[d]=1- P[1].x.b[d];
      }

      // Check P2#P1
      for (d=0;d<D;d++) if(P[2].x.b[d]!=  P[1].x.b[d]) goto global3;
      continue; // They are the same
   global3:

   nn=nn+1;  
      // Compute dist(x,solution)
       distx=distance(1,0,-2);
      // Compute dist(y,solution)
      disty=distance(2,0,-2);
          
      // Compute x fitness
      f[1]=perf(1,function);
      // Compute y fitness
      f[2]=perf(2,function);

      // Check the hypothesis
 
      fx=fabs(f[1]-f[0]); fy=fabs(f[2]-f[0]);

  if((distx-disty)*(fx-fy)>0)
  {
 //printf("\n %f %f %f %f",distx,disty,fx,fy);
    true=true+1;continue;
    }
 // if(distx==disty) true=true+0.5;
  }
   true=0.5*true/nn;

 printf("\n%i %f",radius,true);
 fprintf(f_trace,"%i %i %f %i %f \n", function, D, rate_sample, radius, true);
  }
  break;
  
 case 1: // For local search
 for (radius=2;radius<=D;radius++)
 {
    nn=0;
 for (n=0;n<n_sample;n++)
 {
      for (d=0;d<D;d++) P[0].x.b[d]=alea_integer(0,1); // Random point 0
      // Random point 1
      P[1]=P[0];
      dist=alea_integer(1,radius);
      for (m=0;m<dist;m++)

      {
          d=alea_integer(0,D-1);
          P[1].x.b[d]=1- P[0].x.b[d];
      }
      // Random point 2
      P[2]=P[0];
      dist=alea_integer(1,radius);
      for (m=0;m<dist;m++)
      {

          d=alea_integer(0,D-1);
          P[2].x.b[d]=1- P[0].x.b[d];
      }

      // Check P2#P1
      for (d=0;d<D;d++) if(P[2].x.b[d]!=  P[1].x.b[d]) goto local3;
      continue; // They are the same

      local3:
      nn=nn+1;
      for(m=0;m<3;m++) // Evaluations
        f[m]=perf(m,function);

      distx=distance(1,0,-2);
      disty=distance(2,0,-2);

      fx=fabs(f[1]-f[0]); fy=fabs(f[2]-f[0]);

    if((distx-disty)*(fx-fy)>0) {true=true+1;continue;}
   // if(distx==disty) true=true+0.5;
  }
   true=true/nn;
 printf("\n%i %f",radius,true);
 fprintf(f_trace,"%i %i %f %i %f \n", function, D, rate_sample, radius, true);
  }

 break;
 case 2:   // For local search. Derivation 11 simulation
  for (radius=2;radius<=D;radius++)

 {
    nn=0;
 for (n=0;n<n_sample;n++)
 {
      for (d=0;d<D;d++) P[0].x.b[d]=alea_integer(0,1); // Random point 0
      // Random point 1
      P[1]=P[0];
    //  dist=alea_integer(1,radius);
    //  for (m=0;m<dist;m++)

    rho=alea(1,radius);
    for(m=0;m<rho;m++)
      
      {
          d=alea_integer(0,D-1);
          P[1].x.b[d]=1- P[0].x.b[d];
      }

      for (m=0;m<2;m++) f[m]=perf(m,function); // Evaluations
       s0=0;s1=1;
       if (f[1]>f[0]) {s0=1;s1=0;}   // s1 is the best  (seen as "g")
  
      // Random point 2 around the best
      P[2]=P[s1];
      dist=alea_integer(1,radius);
      for (m=0;m<dist;m++)
      {
          d=alea_integer(0,D-1);
          P[2].x.b[d]=1- P[s1].x.b[d];
      }

      nn=nn+1;

    f[2]=perf(2,function);
    // If the new point is better ...
    if(f[2]<f[s0]) {true=true+1; continue;}
  //  if(f[2]==f[s0]) {true=true+0.5; continue;}
  }
   true=true/nn;
 printf("\n%i %f",radius,true);
 fprintf(f_trace,"%i %i %f %i %f \n", function, D, rate_sample, radius, true);
  }

 break;

 case 3:  // Local search just around the current position
   for (radius=2;radius<=D;radius++)
 {
    nn=0;
 for (n=0;n<n_sample;n++)
 {
      for (d=0;d<D;d++) P[0].x.b[d]=alea_integer(0,1); // Random point 0
      // Random point 1
      P[1]=P[0];
      dist=alea_integer(1,radius);
      for (m=0;m<dist;m++)
      {
          d=alea_integer(0,D-1);
          P[1].x.b[d]=1- P[0].x.b[d];
      }

      for (m=0;m<2;m++) f[m]=perf(m,function); // Evaluations
 
      nn=nn+1;

      
    // If the new point is better ...
    if(f[1]<f[0]) {true=true+1; continue;}
    //if(f[2]==f[s0]) {true=true+0.5; continue;}
  }
   true=true/nn;
 printf("\n%i %f",radius,true);
 fprintf(f_trace,"%i %i %f %i %f \n", function, D, rate_sample, radius, true);
  }

 break;
 } // end switch option
 
  end:;
}
//===========================================================
void init_swarm(int option)
{  /*
 option 0 => positions completely at random
 option 1 => (removed. Was not good)
 option <0 => read from file init.txt
*/
 double CDk[S_max];
int check;
int d;
double dd;
int dist0[S_max],dist1[S_max];
int Dist[S_max][S_max];
//int dmin[1000];// Just for information
int H[D_max]; // Just for information about distances
int k,kk,kkk;
int k1,k2;
double mean0,mean1;
double max;
int min;
double min_r;
int Nk[D_max];
//int norm;
double r;
double R[S_max][D_max];
int St;
int s,s1,s2;
struct position sample[2000];
int sample_s[2000];
int sample_n[S_max];
double sample_d[2000];
double sample_R[2000][D_max];
int space;
double tot;
 double var0,var1;
 
 for (s=0;s<S;s++) P[s].x.size=D;

switch(option)
{
  default: // Completely at random
for (s=0;s<S;s++)
{
	for (d=0;d<D;d++) P[s].x.b[d]=alea_integer(0,1);
}
break;
 


case 2: // Regular distribution

for(d=0;d<D;d++) P[0].x.b[d]=0;  // First position 000...000
//for(d=0;d<D;d++) P[1].x.b[d]=1;  // Second position 111...111
St=0;
loop2:
St=St+1;   // New position rank
for (d=0;d<D;d++)
{
   // Try 0
	P[St].x.b[d]=0;
  for (s=0;s<St;s++)
  {
    dist0[s]=0;
	for (k=0;k<=d;k++)
	{
	  if(P[St].x.b[k] != P[s].x.b[k]) dist0[s]=dist0[s]+1;
	}
   }
   // Mean
    mean0=0;
     for (s=0;s<St;s++) mean0=mean0+dist0[s];
    mean0=mean0/St;
   
   // Variance
    var0=0;
      for (s=0;s<St;s++) var0=var0+ (dist0[s]-mean0)*(dist0[s]-mean0);
   
   // Try 1
	P[St].x.b[d]=1;
   for (s=0;s<St;s++)
  {
 	dist1[s]=0;
	for (k=0;k<=d;k++)
	{
	  if(P[St].x.b[k]!= P[s].x.b[k]) dist1[s]=dist1[s]+1;
	}
   }

   // Mean
   mean1=0;
     for (s=0;s<St;s++) mean1=mean1+dist1[s];
    mean1=mean1/St;
   // Variance
   var1=0;
      for (s=0;s<St;s++) var1=var1+ (dist1[s]-mean1)*(dist1[s]-mean1);


   // Choose the best
	  if(mean1>mean0) {P[St].x.b[d]=1; continue;}
	  if(mean1<mean0) {P[St].x.b[d]=0; continue;}

		if(var1<var0) P[St].x.b[d]=1;
		else P[St].x.b[d]=0;

//printf("\nSt %i, %f %f, %i",St,var0,var1,P[St].x.b[d]);
//scanf("%i",&k);

} // end for(d=0


//printf("\n ");for (s=0;s<St;s++) printf(" %i",dist[s]);

 if(St<S-1) goto loop2;
 
//break;
 // Translate the generated positions
// according to a random "vector"
for (d=0;d<D;d++)
{
   k=alea_integer(0,1);
   if(k==0) continue;
	for(s=0;s<S;s++)
	{
		 P[s].x.b[d]=P[s].x.b[d]+k;
		if(P[s].x.b[d]>1) P[s].x.b[d]=0;
	}
}

 break;


 case -1: // read from file

f_init=fopen("f_init.txt","r");
for (s=0;s<S;s++)
{
printf("\n");
	for (d=0;d<D;d++) 
	{
		fscanf(f_init,"%i",&k);
		P[s].x.b[d]=k;
	}

}
fclose(f_init);
break;

}// end switch(option)

// Velocities 
for (s=0;s<S;s++)
{
	V[s].size=D;
	for (d=0;d<D;d++)
    V[s].b[d]=alea_integer(-1,1);
}

/*
// Display positions
printf("\n Initial %i positions",S);
for (s=0;s<S;s++)
{
  printf("\n");
  for (d=0;d<D;d++) printf("%i",P[s].x.b[d]);
}
*/

  // Save positions
fprintf(f_trace,"\n Initial %i positions",S);
for (s=0;s<S;s++)
{
  fprintf(f_trace,"\n");
  for (d=0;d<D;d++) fprintf(f_trace,"%i",P[s].x.b[d]);
}


// Display distance statistics
// (number of each distance value)
for(d=0;d<=D;d++) H[d]=0;
for(s1=0;s1<S;s1++)
for (s2=s1+1;s2<S;s2++)
{
   d=(int)distance(s1,s2,-2);
   H[d]=H[d]+1;
}
printf("\n Distance statistics: ");
for(d=0;d<=D;d++) printf(" %i",H[d]);

 // Check how evenly is the repartition
 // Note that this just an estimation by sampling

 max=0;
 check=100*D; // Arbitrary ....
 
 for(s1=0;s1<check;s1++)
 {
   for (d=0;d<D;d++) P[S].x.b[d]=alea_integer(0,1); // Draw a point at random
   // Find the nearest pivot (particle)
   min=D+1;
   for(s2=0;s2<S;s2++)
   {
        d=(int)distance(s2,S,-2);
        if(d<min) min=d;
   }
   if(min>max) max=min;
 }

 printf("\nAfter init");
 printf("\nMax radius %i",(int)max);
 //scanf("%i",&d);
}
//===========================================================
double number(int d,struct bits x,double scale)

{
// Compute the value (base 10) of the bit string

  int D;
  double z;
  D=x.size;
  if (d==D-1) return x.b[D-1]*scale;
  z=x.b[d]*scale+2*number(d+1,x,scale);

  return z;
}
//===========================================================
double perf(int s,int function)
{/* Evaluate the function value at position x
 For some functions, there is a trick to save some computational time:
 when the function value is progressively computed as a sum of non negative values
 there is usually no need to compute it entirely. You just have to be sure it does
 not improve the best solution known by the particle

*/
int	bit_sum;
struct bits bits;
int DD;
int    d;
int		dist;
double  f,f1,f2;
int		i,j,k;
double max;
int   num;
double scale;
struct bits x;

x=P[s].x;


nb_eval=nb_eval+1;

switch (function)
{

case 1: // Goldberg's order-3problem:
	// Be sure dimension=3*k
	// The min value is 0  on 1111...111
f=0;
for (d=0;d<D;d=d+3)
{
	bit_sum=0;
	for (i=d;i<d+3;i++) bit_sum+=x.b[i];

  if (bit_sum==0) {f+=0.9; continue;}
 //  if (bit_sum==0) {f+=0.1; continue;} 
	if (bit_sum==1) {f+=0.6; continue;}
	if (bit_sum==2) {f+=0.3; continue;}
	if (bit_sum==3) {f+=1.0; continue;}
}

f=D/3-f;
break;


case 2: // Bipolar order-6
        // Min 0 on any position ..6x0..6x1..
 f=0;

for (d=0;d<D-3;d=d+6)
{
	bit_sum=0;
	for (i=d;i<d+6;i++) bit_sum+=x.b[i];

	if (bit_sum==0) {f+=1.0; continue;}
	if (bit_sum==1) {f+=0.0; continue;}
	if (bit_sum==2) {f+=0.4; continue;}

	if (bit_sum==3) {f+=0.8; continue;}
  if (bit_sum==4) {f+=0.4; continue;}
  if (bit_sum==5) {f+=0.0; continue;}
  if (bit_sum==6) {f+=1.0; continue;}
}

f=D/6-f;
 break;


case 3: // Muhlenbein's order 5
        // Min 0 on 000000000...
max=-P_m[s].f.fi[0]+4 ;

f=0;
bits.size=5;
for (d=0;d<D-2;d=d+5)
{
	for (i=d;i<d+5;i++) bits.b[i-d]=x.b[i];
   num=(int)number(0,bits,1);

   if (num==0) {f+=4.0; goto d3;}
   if (num==1) {f+=3.0; goto d3;}
   if (num==3) {f+=2.0; goto d3;}
   if (num==7) {f+=1.0; goto d3;}
   if (num==31) {f+=3.5; goto d3;}
d3:

  if (init==1) continue;

   if (f<=max+(double)(4*d)/5) // No hope to do better
   {
     save_time=save_time + (double)(D-d-5)/D;
    goto end3;
   } 
}
end3:

f=4*D/5-f;
break;

case 4: // Clerc's order 3 problem 1
	// The min value is 0 on 011*
f=0;
for (d=0;d<D;d=d+3)
{
	bit_sum=x.b[d+2];
	for (i=d;i<d+2;i++) bit_sum+=(int)(x.b[i]*pow(2,2+d-i));
  if (bit_sum==0) {f+=0.9; continue;} //0 =000
	if (bit_sum==1) {f+=0.6; continue;} //1 =001
	if (bit_sum==2) {f+=0.3; continue;} //2 =010
	if (bit_sum==3) {f+=1.0; continue;} //3 =011

 	if (bit_sum==4) {f+=0.2; continue;} //4 =100

  if (bit_sum==5) {f+=0.4; continue;} //5 =101
  if (bit_sum==6) {f+=0.6; continue;} //6 =110
  if (bit_sum==7) {f+=0.8; continue;} //7 =111

}
f=D/3-f;
break;


case 5: // Clerc's order 3 problem 2
	// The min value is 0 on 001*
f=0;
for (d=0;d<D;d=d+3)
{
	bit_sum=x.b[d+2];
	for (i=d;i<d+2;i++) bit_sum+=(int)(x.b[i]*pow(2,2+d-i));
  if (bit_sum==0) {f+=0.2; continue;} //0 =000
	if (bit_sum==1) {f+=1.0; continue;} //1 =001

	if (bit_sum==2) {f+=0.3; continue;} //2 =010
	if (bit_sum==3) {f+=0.8; continue;} //3 =011
 	if (bit_sum==4) {f+=0.6; continue;} //4 =100

  if (bit_sum==5) {f+=0.4; continue;} //5 =101
  if (bit_sum==6) {f+=0.6; continue;} //6 =110
  if (bit_sum==7) {f+=0.9; continue;} //7 =111

}
f=D/3-f;
break;

case 6: // Clerc's order-3 problem 3 (Zebra3)
  // Min 0 on 000111000111...
max=-P_m[s].f.fi[0]+1 ;
f=0;
num=0;

for (d=0;d<D;d=d+3)
{
  bit_sum=0;
  num=num+1;
	for (i=d;i<d+3;i++) bit_sum+=x.b[i];

  if(num%2==0)
  {
    if (bit_sum==0) {f+=0.9; goto d6;}
	  if (bit_sum==1) {f+=0.6; goto d6;}
	  if (bit_sum==2) {f+=0.3; goto d6;}
	  if (bit_sum==3) {f+=1.0; goto d6;}
   }
   else
   {
    if (bit_sum==0) {f+=1.0; goto d6;}
	  if (bit_sum==1) {f+=0.3; goto d6;}
	  if (bit_sum==2) {f+=0.6; goto d6;}
	  if (bit_sum==3) {f+=0.9; goto d6;}
   }
d6:
if (init==1) continue;
 
   if (f<=max+d/3) // No hope to do better
   {
     save_time=save_time + (double)(D-d-3)/D;
//if (d>15) printf("\n%i %i",n_iter,d);
    goto end6;
   } 
}
end6:
f=D/3-f;

break;

case 7: // Parabola
// Warning. Max D is depending on the machine. Typically 50 for a 32 bits one
// Min 0 on pow(2,D-1)= bit sequence 1<D-2 times 0>
scale=1;
f1=pow(2,D-1)*scale;
f2=number(0,x,scale);
f=f2-f1;
 f=f*f;
break;


case 8: // Using Multimodal generated landscape
       f = 0.0;
       for (j = 0; j < peaks; j++)
       {
	        f1 = 0.0;
	        for (k = 0; k < D; k++)
          {
	           if ((int)x.b[k] == landscape[j][k]) f1++;
//printf("\n%i %i %f", (int)x.b[k],landscape[j][k],f1);
	        }


	        if (f1 > f) f = f1; 

       }

//end8:
       f=1-f/(double)D;
break;

case 9: // Clerc's Trap
  // Min 0 at 0000...000
    bit_sum=0;
  for (d=0;d<D;d++)   bit_sum+=x.b[d];
   if(bit_sum==0) f=0;
   else
    f=1.0/bit_sum;
  break;

case 10: // Whitney DF2

// Min 0 at 111...1111
 // Best parameters seem to be:
 // method 100, S30 K3, random links, 10000 eval => mean eps=3.06
 max=P_m[s].f.fi[0] ;
f=0;
bits.size=4;
for (d=0;d<D;d=d+4)
{
	for (i=d;i<d+4;i++) bits.b[i-d]=x.b[i];
   num=(int)number(0,bits,1);

   if (num==15) {f+=0.0; goto d10;}
   if (num==0) {f+=2.0; goto d10;}
   if (num==1) {f+=4.0; goto d10;}

   if (num==2) {f+=6.0; goto d10;}
   if (num==4) {f+=8.0; goto d10;}
   if (num==8) {f+=10.0; goto d10;}
   if (num==3) {f+=12.0; goto d10;}
   if (num==5) {f+=14.0; goto d10;}
   if (num==6) {f+=16.0; goto d10;}
   if (num==9) {f+=18.0; goto d10;}
   if (num==10) {f+=20.0; goto d10;}
   if (num==12) {f+=22.0; goto d10;}

   if (num==14) {f+=24.0; goto d10;}
   if (num==13) {f+=26.0; goto d10;}
   if (num==11) {f+=28.0; goto d10;}
   if (num==7) {f+=30.0; goto d10;}
   
d10:;

  if (init==1) continue;

   if (f>=max) // No hope to do better
   {
     save_time=save_time + (double)(D-d-4)/D;
    goto end10;
   }

}
end10:;
break;

case 11: // Try to optimise a distribution of points
		// in the search space
		// The purpose is to find a good initialisation
		// for other problems
		// There will be Sp particles in the swarm
		// and the dimension of the problem will be DD
		// The _current_ dimension is then Sp*DD
	DD=D/Sp;
f=0;
for (i=0;i<Sp-1;i++) // "particle" i
{
	for (j=i+1;j<Sp;j++)  // "particle j
	{
		// Compute the distance
		dist=0;
		for (d=0;d<DD;d++)
		{
			if (x.b[d+DD*i]!=x.b[d+DD*j]) dist=dist+1;
		}

		// Sum it
			f=f+dist;
	}
}
f=Sp*(Sp-1)*DD-f; // To minimise
return f;
break;

case -1: // Quadratic problem
		// Q has been read on a file
f=0;
for (i=0;i<D;i++)
for (d=0;d<D;d++)
{
  f=f+Q[i][d]*x.b[i]*x.b[d];
}

f=fabs(D+f); // You may have to change it
break;

case 99: // Test

f=0;
for (d=0;d<D;d++) f=f+x.b[d];
   break;

}
//printf("\n f %f",f);
return f;
}

//===========================================================
void swarm_size()
{
 /*
 Completely independent program
 Just to plot some curves DELTA(S) vs S

*/
  int D,d;
  int DELTA[100];
  int S;
  int s_max=60;
  double two_D_S;

  double s_combin[200];
   D=100; // Choose the dimension >=10, please

   // Fill in sum of combin
   s_combin[0]=1;
   for (d=1;d<=D;d++) s_combin[d]=s_combin[d-1]+combin(D,d);

  // for (d=1;d<=D;d++) printf("\n %.0f",s_combin[d]);
   
   // For each swarm size S, find the radius DELTA(S)
   for(S=2;S<=s_max;S++)
   {
     two_D_S=pow(2,D)/S;
   
     for(d=0;d<D;d++)
     {
       if(s_combin[d+1]<two_D_S) continue;
     //DELTA[S]=d+(s_combin[d+1]-two_D_S)/(s_combin[d+1]-s_combin[d]);
//printf("\n%i %f",S, d+ (s_combin[d+1]- two_D_S)/(s_combin[d+1]-s_combin[d]));
      DELTA[S]=(int)(d+0.5);
      goto next_S;
     }
    next_S:;
   }

  for (S=2;S<=s_max;S++) printf("\n %i %i",S,DELTA[S]);
  for (S=2;S<=s_max;S++) fprintf(f_trace,"\n %i %i",S,DELTA[S]);
}

//===========================================================
double	total_error(struct f f)
{
	double	err=0;
	int		i;
	double	z;


	for (i=0;i<f.size;i++)
	{
		z=f.fi[i]; err=err+z*z;
	}
	return sqrt(err);
}
//===========================================================
//Deterministic pseudo random numbers generator KISS.
/*
 the idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
        x(n)=a*x(n-1)+1 mod 2^32
        y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
        z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127
*/


static ulong kiss_x = 1;
static ulong kiss_y = 2;
static ulong kiss_z = 4;
static ulong kiss_w = 8;

static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;


void seed_rand_kiss(ulong seed)
{
    kiss_x = seed | 1;

    kiss_y = seed | 2;
    kiss_z = seed | 4;

    kiss_w = seed | 8;
    kiss_carry = 0;
}

ulong rand_kiss()
{
    kiss_x = kiss_x * 69069 + 1;
    kiss_y ^= kiss_y << 13;


    kiss_y ^= kiss_y >> 17;
    kiss_y ^= kiss_y << 5;
    kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
    kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
    kiss_z = kiss_w;
    kiss_w = kiss_m;
    kiss_carry = kiss_k >> 30;
    return kiss_x + kiss_y + kiss_w;
}


