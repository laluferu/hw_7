n_steps = 15000 
n_burn = 2000

all: poblaciones.py

poblaciones.pdf: poblaciones.dat 
	python poblaciones.py

poblaciones.dat: mcmc_lotkavolterra.x
	./mcmc_lotkavolterra.x $(n_steps) $(n_burn) > poblaciones.dat

mcmc_lotkavolterra.x: mcmc_lotkavolterra.c
	cc mcmc_lotkavolterra.c -lgsl -lgslcblas -lm -o mcmc_lotkavolterra.x

clean:
	rm -f mcmc_lotkavolterra.x  *.pdf poblaciones.dat 
