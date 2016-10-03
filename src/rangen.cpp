
#include "rangen.h"
#include <sys/time.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *cm_rng;


RanGen::RanGen(){
	struct timeval a;
	struct timezone b;
	unsigned long int seed;
	
	gettimeofday(&a,&b);

    FILE*filein = fopen("seed.in","r");
    FILE*fileout = fopen("seed.out","w");

    if (filein == NULL)   
	{
		printf("#seed.in not found. Creating a new seed...\n");
		seed = a.tv_usec;
		//printf("seed = %ld \n",seed);
		//srandom(seed);
		cm_rng = gsl_rng_alloc(gsl_rng_mt19937); //Makoto
		gsl_rng_set(cm_rng,seed);
		fprintf(fileout,"%ld",seed);
		fclose(fileout);
    }
    else
    {
		printf("\n#reading seed.in .\n");
		fscanf(filein,"%ldu",&seed);
		fprintf(fileout,"%ld",seed);
		//srandom(seed);
		cm_rng = gsl_rng_alloc(gsl_rng_mt19937); //Makoto
		gsl_rng_set(cm_rng,seed);
		fclose(filein);
		fclose(fileout);
    }
}
RanGen::~RanGen(){};

bool RanGen::ranbool(){
	int result= (int) gsl_rng_uniform_int (cm_rng,2);
	return (result == 0)? false : true;
}

int RanGen::ranval(int min, int max){
	int diff=max-min;
	int offset= (int) gsl_rng_uniform_int (cm_rng,diff);
	return min + offset;
}

double RanGen::randouble(){
	return gsl_rng_uniform_pos (cm_rng);
}


double RanGen::rangauss(double sigma){
	return gsl_ran_gaussian( cm_rng, sigma);
}


int RanGen::ranpoisson(double mu){
	return gsl_ran_poisson( cm_rng, mu);
}
		
double RanGen::ranpoissonpdf(unsigned int k, double mu){
	return gsl_ran_poisson_pdf( k, mu);
}
