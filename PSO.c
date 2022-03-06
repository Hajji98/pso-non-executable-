struct result PSO (struct param param, struct problem pb, int indFunc) 
{  
	int coop, coop1, coop2;
	int d; 
	struct intDouble intDouble[S_max];
	double error;   
	double errorPrev;
	int followInit=follow1;// For information for further analysis
	// of the variable topology
	int followStop=0; 

	int g;  
	struct position Gr;
	struct vector GX={0}; 

	double helpMin;
	int helpUpdate;
	int index[S_max];     
	int iter; 		// Iteration number (time step)
	int m; 

	//---For TEST. Hard coded parameters for mutation (if option b used in mutate() )
	double minMut=0.7; //0.75;
	double maxMut=1.3; //1.25;
	//---

	int n;
	int noStop;
	double out;
	double p;
	double probaHelp[S_max];
	struct vector PX;	
	struct result R;
	double rad;
	int removeNb; // In case of adaptation
	int s0, s,s1, sd; 
	struct vector V1,V2; 
	double w1,w2,w3;
	double xMax,xMin;
	double zz;

	PX.size=pb.SS.D; 
	GX.size=pb.SS.D; 
	V1.size=pb.SS.D;
	V2.size=pb.SS.D;
	Gr.size=pb.SS.D;
	// -----------------------------------------------------
	// INITIALISATION
	p=param.p; // Probability threshold for random topology
	R.SW.S = param.S; // Size of the current swarm
	removeNb=param.adaptSelect*param.S; // Number of particles to replace
	// (in case of adaptation)

	// Coefficients. 
	for(s=0;s<R.SW.S;s++)
	{ R.SW.w[s]=param.w;
		R.SW.c[s][0]=param.c1;
		R.SW.c[s][1]=param.c2;
	}

	// Help matrix (who has helped who, and when)
	for (s0 = 0; s0 < R.SW.S; s0++)
	{
		for (s=0;s<R.SW.S;s++)
		{
			R.SW.help[s0][s]=param.helpInit; 		
		}
	}

	// Reputations
	for (s0 = 0; s0 < R.SW.S; s0++) R.SW.reput[s0]=param.reputInit;

	// Position and velocity
	for (s = 0; s < R.SW.S; s++)   
	{
		R.SW.X[s].size = pb.SS.D;
		R.SW.V[s].size = pb.SS.D;
	}

	if(pb.SS.normalise>0) {xMin=0; xMax=pb.SS.normalise;} // [0,normalise]^D


	for (s = 0; s < R.SW.S; s++)   
	{

		for (d = 0; d < pb.SS.D; d++)  
		{  
			if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}
			R.SW.X[s].x[d] = alea (xMin,xMax,param.randCase);		
			R.SW.V[s].v[d] = alea( xMin-R.SW.X[s].x[d], xMax-R.SW.X[s].x[d],param.randCase ); 
			// So that  xMin<= x+v <= xMax
		}
	}

	// Take quantisation into account
	if(pb.SS.quantisation==1)
	{
		for (s = 0; s < R.SW.S; s++)
		{  
			R.SW.X[s] = quantis (R.SW.X[s], pb.SS);
		}
	}

	// First evaluations
	errMax=0;
	errMin=infinity;
	for (s = 0; s < R.SW.S; s++) 
	{	
		R.SW.X[s].f = perf (R.SW.X[s], pb.function,pb.SS,pb.objective);
		R.SW.P[s] = R.SW.X[s];	// Best position = current one
	}

	//
	// If the number max of evaluations is smaller than 
	// the swarm size, just keep evalMax particles, and finish
	if (R.SW.S>pb.evalMax) R.SW.S=pb.evalMax;	
	R.nEval = R.SW.S;

	// Find the best
	R.SW.best = 0;
	errorPrev =R.SW.P[R.SW.best].f; // "distance" to the wanted f value (objective)

	for (s = 1; s < R.SW.S; s++)     
	{
		zz=R.SW.P[s].f;
		if (zz < errorPrev)
		{
			R.SW.best = s;
			errorPrev=zz;
		} 
	}

	// Display the best
	//	printf( "\nBest value after init. %1.20e ", errorPrev ); printf(" ");
	//		printf( "\n Position :\n" );
	//		for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );

	// ---------------------------------------------- ITERATIONS	
	noStop = 0;	
	error=errorPrev;	
	iter=0; 
	helpUpdate=0;
	if(param.coop<=5) coop=param.coop;
	else
	{
		if(param.coop<99) // Hybrid nm
		{
			coop1=param.coop/10;
			coop2=param.coop-10*coop1;
			coop=coop1;
		}
		else coop=1; // Cyclically 1 to 4 if no improvement
	}

	while (noStop == 0) 
	{	
		/*
		 // Display the swarm		
		 printf("\n Positions (%i) \ Velocities (%i) after iter %i.",pb.SS.D,pb.SS.D, iter-1 );
		 for (s = 0; s < R.SW.S; s++) 
		 {
			 printf("\n");
			 for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.X[s].x[d] );
			 printf("\ ");
			 for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.V[s].v[d] );
	}
	*/

		// Save the adaptation, for further analysis
		if(coop==3)
		{
			fprintf(f_trace,"%i ", iter );
			for(s=0;s<R.SW.S;s++)
			{			
				fprintf(f_trace,"%f ",R.SW.w[s]);
			}
			fprintf(f_trace,"%i ", iter );
			for(s=0;s<R.SW.S;s++)
			{			
				fprintf(f_trace,"%f ",R.SW.c[s][0]);
			}
			fprintf(f_trace,"%i ", iter );
			for(s=0;s<R.SW.S;s++)
			{			
				fprintf(f_trace,"%f ",R.SW.c[s][1]);
			}
			fprintf(f_trace,"\n");
		}

		switch(param.randNum)
		{
			default:
			case 0:
				for (s = 0; s < R.SW.S; s++)  index[s]=s;  // No random permutation
				break;

			case 1:
				aleaIndex(index, R.SW.S,param.randCase); // Random numbering of the particles
				break;
		}	

		iter=iter+1;
		// Loop on particles
		for (s0 = 0; s0 < R.SW.S; s0++)	// For each particle ...
		{		
			s=index[s0];

			if(iter>1 && (param.helpUpdate==0 || (param.helpUpdate==1 && helpUpdate==1)))
			{
				switch(coop) // Help
				{
					case 1: // Probabilistic reciprocity
						// The more times, and the more recently s has helped g, 
						// the more g will tend to help s
						switch(param.helpOption)
					{
						case 0: // norm a
							zz=0;
						for(g=0;g<R.SW.S;g++) zz=zz+R.SW.help[g][s];

						if(zz>0)
						{
							for(g=0;g<R.SW.S;g++) 
							{
								if(R.SW.help[g][s]>param.helpInit) probaHelp[g]=1;
								else
								{
									probaHelp[g]=R.SW.help[g][s]/zz;
									if(probaHelp[g]<param.helpMin) probaHelp[g]=param.helpMin;
								}
							}
						}
						else for(g=0;g<R.SW.S;g++) probaHelp[g]=param.helpMin;
						break;

						case 1: //norm b
							for(g=0;g<R.SW.S;g++) 
						{
							if(R.SW.help[g][s]>param.helpInit) probaHelp[g]=1;
							else probaHelp[g]=param.helpMin; 
						}	
						break;
						
						case 2: //norm c. Adaptive helpMin
							// Count the number of sure informants (total = sd)
							sd=0;
						for(g=0;g<R.SW.S;g++)  
						{
							if(R.SW.help[g][s]>param.helpInit) 
							{probaHelp[g]=1; sd++;}
						}	
						// If sd<D, assign non null helpMin to the S-sd particles
						helpMin=hminAdapt(param.alpha, R.SW.S-sd, helpSubOption);

						for(g=0;g<R.SW.S;g++) 
						{
							if(R.SW.help[g][s]<=param.helpInit) 
								probaHelp[g]=helpMin;
						}			

						break;

						case 10: // In this case, the help is already a probability
							for(g=0;g<R.SW.S;g++) probaHelp[g]=R.SW.help[g][s];
						break;
					}

					case 2: // Neighbours (Vicinity)
						g=0;
					for(m=0;m<R.SW.S;m++) // All distances
					{
						if(m==s) continue;
						intDouble[g].d=distanceL(R.SW.X[s],R.SW.X[m],2);
						intDouble[g].i=m;
						g++;
					}
					// Sort
					qsort(intDouble,R.SW.S-1,sizeof(intDouble[0]), compareIntDouble); 

			
					// Assign non null probability to the D first ones
					d=pb.SS.D; if (d>R.SW.S-1) d=R.SW.S-1;
					for(m=0;m<d;m++) probaHelp[intDouble[m].i]=1;

					// Assign smaller probability to the other ones
					if(param.helpOption==2)
						helpMin=hminAdapt(param.alpha, R.SW.S-d, helpSubOption);
					else helpMin=param.helpMin;
	
					if(d<R.SW.S-1)
					{
						for(m=d;m<R.SW.S-1;m++) probaHelp[intDouble[m].i]=helpMin; // b		
					}

					break;

					case 3: // Genetically similar (Kin)
						// WARNING. Needs mutation-selection process
						// Compute the "genetical" distance, and sort
						g=0;
					for(m=0;m<R.SW.S;m++) // All distances
					{
						if(m==s) continue;
						intDouble[g].d=distanceGen(R.SW.w[s],R.SW.c[s][0], R.SW.c[s][1],
						                           R.SW.w[m],R.SW.c[m][0], R.SW.c[m][1],
						                           param.w,param.c1,param.c2);
						intDouble[g].i=m;
						g++;
					}
					// Sort
					qsort(intDouble,R.SW.S-1,sizeof(intDouble[0]), compareIntDouble); 

					// Assign non null probability to the D first ones
					d=pb.SS.D; if (d>R.SW.S-1) d=R.SW.S-1;
					for(m=0;m<d;m++) probaHelp[intDouble[m].i]=1; 

					// Assign smaller probability to the other ones
					if(d<R.SW.S-1)
					{
						for(m=d;m<R.SW.S-1;m++) probaHelp[intDouble[m].i]=0; //aa
					}

					break;

					case 4: // Reputation
						// Normalise. Transform the Reputation into probability
						switch(param.helpOption)
					{
						// Normalise. Transform the Reputation into probability
						case 0: // norm a
							zz=0;
						for(g=0;g<R.SW.S;g++) zz=zz+R.SW.reput[g];

						if(zz>0)
						{
							for(g=0;g<R.SW.S;g++)
							{
								if(R.SW.reput[g]>param.reputInit) probaHelp[g]=1;
								else
								{
									probaHelp[g]=R.SW.reput[g]/zz;
									if(probaHelp[g]<param.reputMin) probaHelp[g]=param.reputMin;
								}
							}
						}
						else for(g=0;g<R.SW.S;g++) probaHelp[g]=param.reputMin;
						break;

						case 1: // norm b
							for(g=0;g<R.SW.S;g++) 
						{
							if(R.SW.reput[g]>param.reputInit) probaHelp[g]=1;
							else probaHelp[g]=param.reputMin; 
						}	
						break;

						default:
							case 2: // norm c. Adaptive helpMin
								// Count the number of sure informants (total = sd)
								sd=0;
							for(g=0;g<R.SW.S;g++)  
						{
							if(R.SW.reput[g]>param.reputInit) 
							{probaHelp[g]=1; sd++;}
						}	
							// If sd<D, assign non null helpMin to the S-sd particles
								helpMin=hminAdapt(param.alpha, R.SW.S-sd, helpSubOption);
							for(g=0;g<R.SW.S;g++) 
						{
							if(R.SW.reput[g]<=param.reputInit) 
								probaHelp[g]=helpMin;
						}			

							break;

							case 10: // In this case, the Reputation is already in [0,1]
								for(g=0;g<R.SW.S;g++) probaHelp[g]=R.SW.reput[g];
							break;
					}
					break;
					
					case 5: // Anybody. Any g would accept to help s
						for(g=0;g<R.SW.S;g++) probaHelp[g]=1;
					break;
				}
			}

			// Select the best informant
			g=s;
			for(m=0;m<R.SW.S;m++)
			{
				if(alea(0,1,param.randCase)>probaHelp[m]) continue;
				// m is OK to help s
				if(R.SW.P[m].f<R.SW.P[g].f) g=m; 

				// Save this information for further analysis of the variable topology
				if(m==followInit && m!=s && followStop==0 )
				{
					//fprintf(f_proba,"%i %i %i\n", iter, m,s);
					followInit=s;
					if(followInit==follow2) 
					{
						followStop=1;
						//fprintf(f_proba,"%i %i %i\n", iter, follow1, follow2); // Time to be informed
						if(iter<followMax) followIter[iter][indFunc]++; 
						else 
						{
							printf("\nWarning, iter %i >= followMax %i",iter, followMax);
						}
						//printf("\niter %i", iter);
					}
				}
			}

			// Update the help matrix (g is helping s)
			switch(param.updateOption)
			{
				default: // Linear
					R.SW.help[g][s]=R.SW.help[g][s]+param.helpStep;
				break;

				case 1: // Proportional (logistic)
					R.SW.help[g][s]=R.SW.help[g][s]+param.helpStep*(1-R.SW.help[g][s]);
				break;
			}

			//.. compute the new velocity, and move

			// Exploration tendency
			for (d = 0; d < pb.SS.D; d++)
			{
				R.SW.V[s].v[d]=R.SW.w[s] *R.SW.V[s].v[d];
			}

			// Prepare Exploitation tendency  p-x
			for (d = 0; d < pb.SS.D; d++)
			{
				PX.v[d]= R.SW.P[s].x[d] - R.SW.X[s].x[d];
			}		


			for (d = 0; d < pb.SS.D; d++) 
			{
				GX.v[d]= R.SW.P[g].x[d] - R.SW.X[s].x[d];
			}

			// Gravity centre Gr
			w1=1; w2=1; w3=1; 

			if(g==s) w3=0;

			zz=w1+w2+w3;
			w1=w1/zz; w2=w2/zz; w3=w3/zz;

			for (d = 0; d < pb.SS.D; d++) 
			{
				Gr.x[d]=w1*R.SW.X[s].x[d] + w2*(R.SW.X[s].x[d] + R.SW.c[s][0]*PX.v[d]);
				if(g!=s) Gr.x[d]=Gr.x[d]+ w3*(R.SW.X[s].x[d] +R.SW.c[s][1]*GX.v[d]) ;	

				V1.v[d]= Gr.x[d]-R.SW.X[s].x[d]; // Vector X-G
			}

			// Random point around
			rad=distanceL(R.SW.X[s],Gr,2);

			V2=alea_sphere(pb.SS.D,rad,param.distrib, param.distribMean,param.distribSigma,param.randCase); 

			// New "velocity"
			for (d = 0; d < pb.SS.D; d++) 
			{
				R.SW.V[s].v[d]=R.SW.V[s].v[d]+V1.v[d] + V2.v[d]; 
			}

			// New position
			for (d = 0; d < pb.SS.D; d++) 
			{	
				R.SW.X[s].x[d] = R.SW.X[s].x[d] + R.SW.V[s].v[d];			
			}

			if (R.nEval >= pb.evalMax)  goto end;

			// --------------------------

			if(pb.SS.quantisation==1)
				R.SW.X[s] = quantis (R.SW.X[s], pb.SS);

			// Confinement			
			out=0;
			switch(param.confin)
			{
				default: 
					for (d = 0; d < pb.SS.D; d++)
				{	
					if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}

					if (R.SW.X[s].x[d] < xMin)
					{	
						R.SW.X[s].x[d] = xMin;
						R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
						out=1;
					}
					else
					{
						if (R.SW.X[s].x[d] > xMax)
						{			
							R.SW.X[s].x[d] = xMax;
							R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
							out=1;
						}
					}				
				}	
					break;

					case 1: // No confinement and no evaluation if outside ("let if fly")
						for (d = 0; d < pb.SS.D; d++)
				{	
					if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}

					if (R.SW.X[s].x[d] < xMin || R.SW.X[s].x[d] > xMax) out=1;		
				}	
					break;
			}

			if(pb.SS.quantisation==1 && out>0)
			{
				R.SW.X[s] = quantis (R.SW.X[s], pb.SS);
			}		
			// If the position is inside
			if(param.confin==0 || (param.confin==1 && out<zero))
			{
				// Evaluation
				R.SW.X[s].f =perf(R.SW.X[s],pb.function, pb.SS,pb.objective);
				R.nEval = R.nEval + 1;				
				// ... update the best previous position		
				if (R.SW.X[s].f < R.SW.P[s].f)	// Improvement
				{															
					R.SW.P[s] = R.SW.X[s]; // Replace the previous best

					switch(param.updateOption) // Increase reputation
					{
						default: 
							case 0: // Linear
							R.SW.reput[s]=R.SW.reput[s]+param.reputStep; 
							break;
						case 1: // Proportional (logistic)
							R.SW.reput[s]=R.SW.reput[s]+param.reputStep*(1-R.SW.reput[s]);
							break;
					}
					// ... update the best of the bests
					if (R.SW.P[s].f < R.SW.P[R.SW.best].f)
					{		
						R.SW.best = s;			
					}			
				}		
			}

			if(param.trace>0)
			{
				// Keep trace of every position, for further tests
				fprintf(f_trace,"%i %f ",s,R.SW.X[s].f);
				for (d = 0; d < pb.SS.D; d++)
				{
					fprintf(f_trace,"%f ",R.SW.X[s].x[d]);
				}
				fprintf(f_trace,"\n");
			}	
		}			// End of "for (s0=0 ...  "	

		// Check if finished
		error = R.SW.P[R.SW.best].f;

		if (error < errorPrev)	// Improvement of the global best
		{		
			helpUpdate = 1;	// Help matrix will not be modified
			// (if param.helpUpdate==1)
		}
		else			// No global improvement
		{			
			helpUpdate = 0;	// Help matrix will be modified	

			if(param.coop>5) // Hybrid
			{
				if(param.coop<99)
				{
					if(coop==coop1) coop=coop2;
					else coop=coop1;
				}
				else
				{
					coop=coop+1; if(coop==5) coop=1; // Cycle 1 to 4
				}
			}
		}

		errorPrev = error;
		end:
			/*
			 // Save reputations
			 for(m=0;m<R.SW.S;m++)
			 fprintf(f_trace,"%.1f ",R.SW.reput[m]);
			 fprintf(f_trace,"\n");
			 */
			if (error > pb.epsilon && R.nEval < pb.evalMax)
		{
			noStop = 0;	// Won't stop

			// Decay help
			for(s0=0;s0<R.SW.S;s0++)
			{
				switch(param.updateOption)
				{ 
					default:
						case 0: // 0. Linear
						for(s1=0;s1<R.SW.S;s1++)
					{
						zz=R.SW.help[s0][s1]-param.helpDecay;
						if(zz<0)	R.SW.help[s0][s1]=0; 
						else R.SW.help[s0][s1]=zz;
					}
					break;
					
					case 1: // Proportional
						for(s1=0;s1<R.SW.S;s1++)
							R.SW.help[s0][s1]=param.helpDecay*R.SW.help[s0][s1];					
					break;
				}
			}

			// Decay reputations
			for(m=0;m<R.SW.S;m++)
			{
				switch(param.updateOption)
				{
					default: 
						case 0: // Linear
						zz=R.SW.reput[m]-param.reputDecay;
					if(zz<0)  R.SW.reput[m]=0; 
					else R.SW.reput[m]=zz;	
					break;

					case 1: // Proportional
						R.SW.reput[m]=param.reputDecay*R.SW.reput[m];
					break;
				}
			}	

			// If adaptation
			if(param.adapt==1 || coop==3)
			{
				// Sort the swarm by increasing fitness
				for(s=0;s<R.SW.S;s++)
				{
					intDouble[s].i=s;
					intDouble[s].d=R.SW.P[s].f;
				}
				qsort(intDouble,R.SW.S,sizeof(intDouble[0]), compareIntDouble);

				// Duplicate and mutate the best ones
				// and replace the worst ones

				for(s0=0;s0<removeNb;s0++)
				{
					s=intDouble[s0].i; // Particle to duplicate-mutate
					m=intDouble[R.SW.S-s0].i; // Particle to remove

					R.SW.w[m]=mutate(R.SW.w[s], param.adaptMutW,param.randCase,minMut*param.w, maxMut*param.w);
					R.SW.c[m][0]=mutate(R.SW.c[s][0],param.adaptMutC,param.randCase,minMut*param.c1,maxMut*param.c1);
					R.SW.c[m][1]=mutate(R.SW.c[s][1],param.adaptMutC,param.randCase,minMut*param.c2,maxMut*param.c2);	
				}

				// Re-initialise the "help" and the reputation of the new particles
				for (s1=0;s1<R.SW.S;s1++) R.SW.help[m][s1]=param.helpInit;
				R.SW.reput[m]=param.reputInit;
			}
		}
			else
		{
			noStop = 1;	// Will stop
		}		
	} // End of "while nostop ...

	// printf( "\n and the winner is ... %i", R.SW.best );			
	R.error = error;
	return R;  
}

// =============================================================================
double hminAdapt(double alpha, int Ssd, int option)
{
	double hmin;
	if(Ssd<=0) return 1;
	switch(option)
	{
		default:
		case 0:
			hmin=1./Ssd; // b
			break;
		case 1:
			hmin=1-pow(1-alpha,1./Ssd); // c
			break;
	}

	return hmin;
}
// =============================================================================
struct position quantis (struct position x, struct SS SS ) 
{     
	/*
	 Quantisation of a position
	 Only values like x+k*q (k integer) are admissible 
		 */ 
	int d;
	double qd;
	struct position quantx;

	quantx = x;     
	for (d = 0; d < x.size; d++)
	{
		qd = SS.q.q[d];	

		if (qd > zero)	// Note that qd can't be < 0
		{           				
			quantx.x[d] = qd * floor (0.5 + x.x[d] / qd);	    
		}	
	}
	return quantx;    
}
