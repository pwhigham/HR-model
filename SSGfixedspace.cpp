/*
Network Model
Hill & Robertson Paper (1966)
Space implicit. 
9nd Sept., 2018
Author: P A Whigham

Updated: 4/2/2020
This is the steady-state population model using weighted edges
for networks.


Updated: 10/10/2018

Changed to have all inputs defined (but still fixed space)
This will allow both Fig 1 and other measures to be done

We currently have defined: "fixed" networks:
 
Ring, SM , SF 
 
Setup that you read in the space - and can just define the folder where the 1 or more
spaces are loaded.  This allows many SM & SF spaces if required, and a single one for 
the ring.

C compiler:  bcc32c <file>
With space it is approx x3 in terms of time, but still acceptable.

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
/* 0 == AB
   1 == aB
   2 == Ab
   3 == ab
*/
// Here we define the parameters of the model that are "fixed"
// Allows us to use main memory vs stack

#define N 64 // Number of individuals

#define SMALLFITNESS 0.000001
#define FALSE 0
#define TRUE 1

static int pop[N][2];  //  Population with 2 gametes - using same mapping as R programme
static int nextpop[N][2];  // Next generation
static float fitness[N];  // Fitness for N individuals
static float cumsum[N];   // Holds calculation for cumulative sum/sum
static int fixedResult[4]; // <fixAA, fixaa, gens>
                           // where fix* = TRUE/FALSE, 
						   // gens is number generations to fixation
/*************************************************************************************/
/*   NETWORK STRUCTURE 
***************************************************************************************/

static float network[N][N];  
					   
/****************************************** VALUES THAT MAY NEED TO BE CHANGED ********
***************************************************************************************/

static int recombine[4][4] = {
			{0,1,2,1},
			{0,1,0,1},
			{0,3,2,3},
			{2,3,2,3}
};

static float fit_table[4][4];


void set_fitness(float s,float t)
{
	
	fit_table[0][0] = 1.0; // AB AB
	fit_table[0][1] = (1-0.5*s); // AB aB
	fit_table[0][2] = (1-0.5*t); // AB Ab       
    fit_table[0][3] = ((1-0.5*s)*(1-0.5*t)); // AB ab
    fit_table[1][0] = (1-0.5*s); // aB AB  
    fit_table[1][1] = (1-s); // aB aB
    fit_table[1][2] = ((1-0.5*s)*(1-0.5*t));  //aB Ab 
    fit_table[1][3] = ((1-s)*(1-0.5*t)); // aB ab
    fit_table[2][0] = (1-0.5*t); //Ab AB                      
    fit_table[2][1] = ((1-0.5*s)*(1-0.5*t)); //Ab aB
    fit_table[2][2] = (1-t); // Ab Ab
    fit_table[2][3] =((1-0.5*s)*(1-t)); // Ab ab 
    fit_table[3][0] =((1-0.5*s)*(1-0.5*t)); // #ab AB
    fit_table[3][1] =((1-s)*(1-0.5*t));  // ab aB
    fit_table[3][2] =((1-0.5*s)*(1-t)); // ab Ab                      
    fit_table[3][3] =((1-s)*(1-t)); // ab ab

	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			fit_table[i][j]+=SMALLFITNESS;
}
/************************** UTILITIES ***************************/

/*****************************************************************
********* trimstr - trim spaces at the end of a string
*****************************************************************/ 
char *trimstr(char *s)
{
    if(s) { /* Don't forget to check for NULL! */
        while(*s && isspace(*s))
            ++s;
        if(*s) {
            register char *p = s;
            while(*p)
                ++p;
            do {
                --p;
            } while((p != s) && isspace(*p));
            *(p + 1) = '\0';
        }
    }
    return(s);
}
/******************************************************************
*********** intToString - convert an integer to the string
******************************************************************/
char* intToString(int x)
{
	int length = snprintf( NULL, 0, "%d",x);
	char* str = (char *) malloc( length + 1 );
	snprintf( str, length + 1, "%d", x );
	return(str);
}
/*
randrange(lower,upper)
Generate random integers in range lower to upper inclusive.
i.e. randrange(0,5)  generates random samples of 0,1,2,3,4,5

*/
int randrange(int lower, int upper) 
{ 
	return (rand() % 
           (upper - lower + 1)) + lower; 
} 


/**************************** NETWORK FUNCTIONS ******************/

  
/*   loadNetwork(fname) Load a network into network array from file
********************************************************************/

void loadNetwork(char* fname)
{
	FILE *net;
    net = fopen(fname, "r");
	
	if (net==NULL)
	{
		printf("FAILED TO OPEN NETWORK: %s\n",fname);
		exit(0);
	}
	for(int i=0;i<N;++i)
		for(int j=0;j<N;++j)
		{
			if (fscanf(net, "%f", &network[i][j])!=1) printf("FAILED\n");
		}
	fclose(net);
}
/*   printNetwork() - print current network to stdoutj
***************************************************************************/
void printNetwork()
{
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
			printf("%f",network[i][j]);
		printf("\n");
	}
}


int mapping(int a1, int a2)
{
  if ((a1==0) && (a2==0)) return(0);
  if ((a1==1) && (a2==0)) return(1);
  if ((a1==0) && (a2==1)) return(2);
  return(3);
}
/*
NOTE: This is only required when setting p or q to zero
Since we don't do this in the model runs, it should make no difference
to the result.  
*/
double rand01()
{
	double res=(((double)rand())/((double)(RAND_MAX)));
	return((res==(double)0.0)?(double)0.0000000001:res); 
}
/*
Initialisation - should be > p rather than >= p ?

*/

void initpop(double p, double q)
{
	int tmp[N][4];
	for (int i=0;i<N;++i)
	{
		tmp[i][0] = (rand01() > p)?1:0;
		tmp[i][2] = (rand01() > p)?1:0;
		tmp[i][1] = (rand01() > q)?1:0;
		tmp[i][3] = (rand01() > q)?1:0;
	}
	for (int i=0;i<N;++i)
	{
		pop[i][0] = mapping(tmp[i][0],tmp[i][1]);
		pop[i][1] = mapping(tmp[i][2],tmp[i][3]);
	}
}
/*

DRIFT ! 

*/

void calc_fitness()
{
//	for (int i=0;i<N;++i)
//		fitness[i] = 1;   /* DRIFT */
	
	
	for (int i=0;i<N;++i)
		fitness[i] = fit_table[pop[i][0]][pop[i][1]];
}
/*
*  fixedAa
	g1 = 1, g2 = 3  is for AA
	g1 = 0, g2 = 2  is for aa
	
*/
int fixedAa(int g1, int g2)
{
	for(int i=0;i<N;++i)
	{
		if ((pop[i][0]==g1) ||
            (pop[i][0]==g2) ||
            (pop[i][1]==g1) ||
			(pop[i][1]==g2))
				return(FALSE);
	}
	return(TRUE);
}
void print_current_pop()
{
	for (int i=0;i<N;++i)
		printf("<IND %d> <G1 %d> <G2 %d>: %f\n",i,pop[i][0],pop[i][1],fitness[i]);
}
/*
Select a parent proportional to fitness from the entire population.
All values need to be calculated since population may change each timestep.
Returns index of parent from population pop[]
Assumes that fitness has been calculated previously.  Note this can be done
once at the beginning after initialisation, then just for the new child (one operation)
See calc_fitness() for function. 
*/
int select_global_parent()
{
	register float total=0.0;
	
	for(register int i=0;i<N;++i)
	{
		total += fitness[i];
		cumsum[i] = total;
	}
	for (register int i=0;i<N;++i)
		cumsum[i]/=total;  // Normalise by the sum of fitness...
	
	/* Now do the prop fitness selection */
	register float rnd = ((float)rand()/(float)(RAND_MAX));
	
	for(register int i=0;i<N;++i)
	{
		if (cumsum[i] >= rnd) return(i);
	}
	printf("ERROR with select_global_parent\n");
	exit(0); /* Should never get here but will crash if it does... */
}


/*
Set the normalised cumulative sum of fitness for location <locn> 
Uses the network[locn][N] - if 1, then include in cumsum running total
*/
void calc_cumsum(int locn)
{
	register float total=0.0;
	
	for(register int i=0;i<N;++i)
	{
		if (network[locn][i] > 0.0)
		{
			// total += fitness[i];   // /* ORIGINAL MODEL */
			total += (fitness[i]*network[locn][i]);  /* Prob. selection fitness * outward edge */
			cumsum[i] = total;
		}
		else cumsum[i] = total;  /* Just set to whatever was before */
	}
	
	for (int i=0;i<N;++i)
		cumsum[i]/=total;  // Normalise by the sum of fitness...
}
/*
Select parent using the cumsum that has been set
*/
int select_parent()
{
	float rnd = ((float)rand()/(float)(RAND_MAX));
	float sum=0.0;
	
	for(int i=0;i<N;++i)
	{
		if (cumsum[i] >= rnd) return(i);
	}
	printf("Failed select parent \n");
	exit(0);
	return(N-1);  // shouldn't ever get here...
}

void calc_locn_cumsum(int locn)
{
	register float total=0.0;
	for(register int i=0;i<N;++i)
	{
		if (network[locn][i] > 0.0)
		{
			total += network[locn][i];  // Assumes equal weighting for all links in graph
			cumsum[i] = total;
		}
		else cumsum[i] = total;  /* Just set to whatever was before */
	}
	for (int i=0;i<N;++i)
		cumsum[i]/=total;  // Normalise by the sum of fitness...
}

/*
IN: r - recombination rate
Create next timestep in steady state model.
Randomly generate the location -- and just do
the one HR step.  
Select both parents from random location deme as per 
generational model.  
Note that last time I found that putting the child at the 
random location didn't work at all.  We will first 
confirm this is the case.
FIRST MODEL:  Try putting child at random location.
RESULTS:  
OUT: Changes just a single member of pop and updates
     the fitness for this location.
	 
*/
void next_timestep(float r)
{
	int p1,p2;  // parents
	int tmp;  
	int parent1[2], parent2[2];
	int childloc;  
	
	
	p1 = select_global_parent();
    calc_cumsum(p1);  /* calculate cumulative fitness in deme for p1 */
	p2 = select_parent();

	parent1[0] = pop[p1][0];
	parent1[1] = pop[p1][1]; // Copy values 
	parent2[0] = pop[p2][0];
	parent2[1] = pop[p2][1];
		
	if (((float)rand()/(float)(RAND_MAX)) < r)
	{
			// recombination of p1
			tmp = recombine[parent1[0]][parent1[1]];
			parent1[1] = recombine[parent1[1]][parent1[0]];
			parent1[0] = tmp;
	}
		
	if (((float)rand()/(float)(RAND_MAX)) < r)
	{
			// recombination of p2
			tmp = recombine[parent2[0]][parent2[1]];
			parent2[1] = recombine[parent2[1]][parent2[0]];
			parent2[0] = tmp;
	}
		// and create the child....
		

	// TRY RANDOMLY SELECTING FROM WITHIN THE DEME...
	
	calc_locn_cumsum(p1);  // Setup for selection from globally selected parent location 
	                       // in the deme of p1
    
	childloc = select_parent();  // and pick the location for child in deme

	pop[childloc][0] = (((float)rand()/(float)(RAND_MAX)) < 0.5)?parent1[0]:parent1[1];
	pop[childloc][1] = (((float)rand()/(float)(RAND_MAX)) < 0.5)?parent2[0]:parent2[1];
	
	// and update the fitness table...
	fitness[childloc] = fit_table[pop[childloc][0]][pop[childloc][1]];

}
/*
 popfixed:  r (recombination)
             p,q - initialisations of population
OP: Run till population fixes on both locus 1 with A or a
    and locus 2 with B or b
OUT: Fills in the static variable fixedResult[4]

*/
void popfixed(float r,float p, float q)
{
	int gens=0;
	int fixA=FALSE, fixB=FALSE;
	int fixAgens, fixBgens;  
	
	initpop((double)p,(double)q);  // p and q 
	calc_fitness();  // calculate fitness for current population
					// Fitness table has been previously created.
					// Just need to do this once because SS looks after single
					// change to this table after creating new individual.

	while(TRUE)
	{
		if (!fixA)
		{
			if (fixedAa(1,3)) // AA - NOTE THIS SEEMS WRONG, but results are correct?
			{
				fixA=1;
				fixAgens=gens;
			}
			else
			if (fixedAa(0,2)) // aa
			{
				fixA=2;
				fixAgens=gens;
			}
		}
		if (!fixB)
		{
			if (fixedAa(2,3)) // BB - NOTE THIS SEEMS WRONG, but results are correct?
			{
				fixB=1;
				fixBgens=gens;
			}
			else
			if (fixedAa(0,1)) // bb
			{
				fixB=2;
				fixBgens=gens;
			}
		}

		if (fixA && fixB) break;
		
		next_timestep(r);
		
		gens+=1;  
	}
	fixedResult[0]=fixA;  
	fixedResult[1]=fixAgens;
	fixedResult[2]=fixB;
	fixedResult[3]=fixBgens; 
}

void panmictic()
{
	for (int i=0;i<N;++i)
		for(int j=0;j<N;++j)
			network[i][j] = 1;
}

void diagzero()
{
	for (int i=0;i<N;++i)
		network[i][i] = 0;
}

int main(int argc, char *argv[])
{
	int numruns,
		numspaces;  // number of spaces in the network folder
					// these are numbered 1 .. numspaces
	int writeheader=1;
	float r,s,t;
	float p,q;
    int Nialpha, Nibeta;
    float Nc;  // converted to r	
	
	char netfname[512]; // pathname for folder holding fixed spaces
	char *networkpath;
		
	
	if (argc<10)
	{
		printf("<p float> <q float> <Nialpha int> <Nibeta int> <Nc> <numruns int> \n <network folder path> <numspaces int> <writeheader int> \n");
		printf("N = %d\n", N);
		return(1);
	}
	sscanf (argv[1],"%f",&p);sscanf (argv[2],"%f",&q);
	sscanf (argv[3],"%d",&Nialpha);sscanf (argv[4],"%d",&Nibeta);
	sscanf(argv[5],"%f",&Nc);
	sscanf(argv[6],"%d",&numruns);
	networkpath = argv[7];
	sscanf (argv[8],"%d",&numspaces);
	sscanf(argv[9],"%d",&writeheader);
	
	srand((unsigned int)time(NULL));	// initialise random number generator
		
	s = (float)Nialpha/(float)N;
	t = (float)Nibeta/(float)N; // See pg 272	
	r = Nc/(float)N;  // set recombination rate
	
	set_fitness(s,t);  // Nialpha, Nibeta calculation - Fixed for Figure 4
	
	if (writeheader) printf("N,Nialpha,Nibeta,Nc,p,q,fixAa,gensAa,fixBb,gensBb\n");

	for (int runs=0;runs < numruns; ++runs)
	{			
		char* strexp = intToString((runs%numspaces)+1);
			
		sprintf(netfname,"%s\\\\%s.txt\n",networkpath,strexp);
		loadNetwork(trimstr(netfname));
		
		//diagzero();
		
		popfixed(r,p,q);
		printf("%d,%d,%d,%f,%f,%f,%d,%d,%d,%d\n",N,(int)Nialpha,(int) Nibeta,
						Nc,p,q,fixedResult[0],fixedResult[1],fixedResult[2],fixedResult[3]);
			
	}
}
