/*
 * mcmc_solar.c
 *
 *  Created on: 15/11/2015
 *      Author: Laura Ferrer
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define USAGE "./mcmc_solar.x n_steps n_burn"
#define PI 3.14159265

int leer(double *year, double *manc, int n_puntos);

int main (int argc, char **argv)
{

	double a = 1;
	double b = 1;
	double c = 1;
	double d = 1;

	double a2;
	double b2;
	double c2;
	double d2;

	int datos = 0;
	int i;
	int iteraciones = atof(argv[2])+atof(argv[1]);
	int  corte = atof(argv[1]);

	srand48(time(NULL));
	FILE *file;
	file = fopen("monthrg.dat", "r");

	while(!feof(file)) //Cuenta la cantidad de datos provistos.
	{
  		char ch = fgetc(file);
  		if(ch == '\n')
  		{
    			datos++;
  		}
	}

	fclose(file);

	double y_real[datos];
	double t[datos];

	datos = leer(t, y_real, datos);

	double y[datos];
	double y2[datos];

	void funcion (int modo)
	{
		if (modo == 1)
		{
			for (i = 0; i < datos; i += 1)
			{
				y[i] =a*cos((2*PI/d)*t[i] + b)+ c;
			}
		}
		if (modo == 2)
		{
			for (i = 0; i < datos; i += 1)
			{
				y2[i] =a2*cos((2*PI/d2)*t[i] + b2)+ c2;
			}
		}
	}

	double parecido (int modo)
	{
		double suma=0;
		if (modo == 1)
		{
			for (i = 0; i < datos; i += 1)
			{	
				suma += pow(y[i]-y_real[i],2.0);
			}
		}
		else if (modo == 2)
		{
			for (i = 0; i < datos; i++)
			{
				suma += pow(y2[i]-y_real[i],2.0);
			}
		}

		return (1.0/2.0)*suma;
	}
	int j;
	for (j = 0; j < datos; j++)
	{	
		if (j>corte){
		srand((unsigned)time(NULL));
		const gsl_rng_type * T;
    		gsl_rng * r;
   		gsl_rng_env_setup();
    		T = gsl_rng_default;
   		r = gsl_rng_alloc (T);

		// cambios al azar
		a2 = a + gsl_ran_gaussian(r,drand48());
		srand((unsigned)time(NULL));
		b2 = b + gsl_ran_gaussian(r,drand48());
		srand((unsigned)time(NULL));
		c2 = c + gsl_ran_gaussian(r,drand48());
		srand((unsigned)time(NULL));
		d2 = d + gsl_ran_gaussian(r,drand48());

		funcion (1);
		funcion(2); //calculo de las nuevas y.


		// Calculo de los errores y similitudes
		double dif1 = parecido(1);
		double dif2 = parecido(2);
		double alpha = dif1/dif2;//Comparación

		if(dif2 < dif1)
		{
			printf ("%g,%g,%g,%g\n", a2, b2,c2, d2);
			a = a2;
			b = b2;
			c = c2;
			d = d2;
		}
		else
		{	srand((unsigned)time(NULL));
			double beta = drand48(); // Otro valor al azar
			if (beta >= alpha)
			{
				printf ("%g,%g,%g,%g\n", a2, b2 ,c2, d2);
				a = a2;
				b = b2;
				c = c2;
				d = d2;
			}
			else
			{
				printf ("%g,%g,%g,%g\n", a, b ,c, d);
			}

		}
	}}
printf("%g %g\n", t[0],y_real[0]);
return 0;
}

int leer(double *year, double *manc, int n_puntos)
{	
	double month[n_puntos];
	
	FILE *file;
	int i;
	file=fopen("monthrg.dat","r");
	rewind(file);
	if(!file)
	{
    		printf("problems opening the file %s\n","monthrg.dat");
    		exit(1);
	}
	double nada1;
	double nada2;

	for(i=0;i<n_puntos;i++)
	{
		do{
			fscanf(file, "%lf %lf %lf %lf %lf\n", &year[i], &month[i], &manc[i], &nada1 , &nada2);

	  	}while(nada1 == -99.0);
	}

	fclose(file);

	int limite = n_puntos-491;
	
	for(i=0;i<limite;i++)
	{
		year[i] = year[i] + (month[i]-1.0)/12.0;
		//printf("%g %g\n", year[i], manc[i]);
	}
	return limite;
}
