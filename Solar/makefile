all: pdf

pdf: dat 
	python solar.py

dat: ejecutable
	./mcmc_solar.x 2000 500 > solar.dat

ejecutable:
	gcc mcmc_solar.c -lgsl -lgslcblas -lm -o mcmc_solar.x

clean:
	rm -f *.pdf *.x solar.dat 
