/*
 * mcmc_lotkavolterra.c
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
#define USAGE "./mcmc_lotkavolterra.x n_steps n_burn"

void leer(double *tiempo, double *presa, double *depredador, int n_puntos);
void iniciar(double *lista, int n_puntos);

int main (int argc, char **argv)
{	
	int i=0;
	int lineas = 0;
	int iteraciones = atof(argv[2])+atof(argv[1]);
	int corte = atof(argv[1]);

	double a = 1;
	double b = 1;
	double c = 1;
	double d = 1;

	double a2;
	double b2;
	double c2;
	double d2;

	FILE *file;
	file = fopen("lotka_volterra_obs.dat", "r");

	while(!feof(file)) //Cuenta la cantidad de datos provistos.
	{
  		char ch = fgetc(file);
  		if(ch == '\n')
  		{
    			lineas++;
  		}
	}

	fclose(file);

	lineas --;// menos la primera línea

	double x[lineas];
	iniciar(x,lineas);
	double y[lineas];
	iniciar(y,lineas);
	double x_real[lineas];
	iniciar(x_real,lineas);
	double y_real[lineas];
	iniciar(y_real,lineas);
	double t[lineas];
	iniciar(t,lineas);
	double lay[lineas];
	iniciar(lay,lineas);
	double lax[lineas];
	iniciar(lax,lineas);

	leer(t,x_real,y_real,lineas);
	lay[0]=y_real[0];
	lax[0]=x_real[0];

	double parecido (double modo)
	{	
		i = 0;
		double suma=0;
		if (modo == 1)
		{

			for (i = 0; i < lineas; i += 1)
			{
				suma += pow(x[i]-y_real[i],2.0);
				suma += pow(y[i]-x_real[i],2.0);
			}

		}

		if (modo == 2)
		{

			for (i = 0; i < lineas; i += 1)
			{
				suma += pow(lax[i]-y_real[i],2.0);
				suma += pow(lay[i]-x_real[i],2.0);
			}
		}

		return (1.0/2.0)*suma;
	}

	double dx (double x, double y, double la, double lb)
	{
		double dx = x*(la-lb*y);
		return dx;
	}

	double dy (double x, double y, double lc, double ld)
	{
		double dy = -y*(lc-ld*x);
		return dy;
	}

	void kutta4 (double la, double lb, double lc, double ld)
	{
		
		double h = t[1]-t[0];
		i = 1;

		for (i = 1; i < lineas; i += 1)
		{
			//Inicio
			double kx1 = dx(lax[i-1],lay[i-1],la,lb);
			double ky1 = dy(lax[i-1],lay[i-1],lc,ld);

			//paso 1
			double x1 = lax[i-1] + (h/2.0)*kx1;
			double y1 = lay[i-1] + (h/2.0)*ky1;

			double kx2 = dx(x1,y1,la,lb);
			double ky2 = dy(x1,y1,lc,ld);

			//paso 2
			double x2 = lax[i-1] + (h/2.0)*kx2;
			double y2 = lay[i-1] + (h/2.0)*ky2;

			double kx3 = dx(x2,y2,la,lb);
			double ky3 = dy(x2,y2,lc,ld);

			//paso 3
			double x3 = lax[i-1] + h*kx3;
			double y3 = lay[i-1] + h*ky3;

			double kx4 = dx(x3,y3,la,lb);
			double ky4 = dy(x3,y3,lc,ld);

			//paso 4

			double mx = (1.0/6.0)*(kx1 + 2.0*kx2+2.0*kx3+kx4);
			double my = (1.0/6.0)*(ky1 + 2.0*ky2+2.0*ky3+ky4);

			lax[i]= lax[i-1] + h*mx;
			lay[i]= lay[i-1] + h*my;

		}
	}

	void mcmc()
	{
		i = 0;
		for (i = 1; i < iteraciones; i++)
		{
			if (i>corte){
    			const gsl_rng_type * T;
    			gsl_rng * r;
   			gsl_rng_env_setup();
    			T = gsl_rng_default;
   			r = gsl_rng_alloc (T);

			a2 = a + gsl_ran_gaussian(r,drand48());
			srand((unsigned)time(NULL));
			b2 = b + gsl_ran_gaussian(r,drand48());
			srand((unsigned)time(NULL));
			c2 = c + gsl_ran_gaussian(r,drand48());
			srand((unsigned)time(NULL));
			d2 = d + gsl_ran_gaussian(r,drand48());
			srand((unsigned)time(NULL));

			kutta(a,b,c,d);
			
			int j =0;

			for (j = 0; j < lineas; j++) // Copia las listas
			{
     				x[j] = lax[j];
				y[j] = lay[j];
   			}

			kutta(a2,b2,c2,d2);

			// Calculo de los errores y similitudes
			double difvieja = parecido(1);
			double difnueva = parecido(2);

			if (difvieja<difnueva)
			{
				a=a2;
        	    		b=b2;
        	   		c=c2;
        	    		d=d2;
				printf ("%g,%g,%g,%g\n",a,b,c,d);
			}
			else
			{	
				srand((unsigned)time(NULL));
				double alfa = difvieja/difnueva;
				double beta = drand48(); 
				
				if (beta >= alfa)
				{
					a=a2;
        	    			b=b2;
        	   			c=c2;
        	    			d=d2;
					printf ("%g,%g,%g,%g\n",a,b,c,d);
				}
				else
				{
					printf ("%g,%g,%g,%g\n",a,b,c,d);
				}	
			}
			}
		}
	}

	return 0;
}

void leer(double *tiempo, double *presa, double *depredador, int n_puntos)
{	
	FILE *file;
	int j;
	file = fopen("lotka_volterra_obs.dat", "r");
	rewind(file);
	for(j=0;j<n_puntos;j++)
	{
		fscanf(file, "%lf %lf %lf\n", &tiempo[j], &presa[j], &depredador[j]);
	}
	fclose(file);
}
void iniciar(double *lista, int n_puntos)
{	
	int k = 0;
	for(k=0;k<n_puntos;k++)
	{
		lista[k] = 0.0;
	}
	
}
