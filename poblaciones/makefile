poblaciones.pdf: poblaciones.dat 
	python poblaciones.py

poblaciones.dat: mcmc_lotkavolterra.x
	./mcmc_lotkavolterra.x > poblaciones.dat

mcmc_lotkavolterra.x: mcmc_lotkavolterra.c
	cc mcmc_lotkavolterra.c -lgsl -lgslcblas -lm -o mcmc_lotkavolterra.x

clean:
	rm -f *.pdf *.dat 
