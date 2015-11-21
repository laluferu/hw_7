all: pdf

pdf: data
	python poblaciones.py

data: ejecutable
	./mcmc_lotkavolterra.x 4000 20 > poblaciones.dat

ejecutable:
	gcc mcmc_lotkavolterra.c -lgsl -lgslcblas -lm -o mcmc_lotkavolterra.x

clean:
	rm -f mcmc_lotkavolterra.x *.pdf poblaciones.dat 
