n_steps = 15000
n_burn = 2000

all: solar.pdf

solar.pdf: solar.py  solar.dat 
	python solar.py

solar.dat: mcmc_solar.x
	./mcmc_solar.x $(n_steps) $(n_burn) > solar.dat

mcmc_solar.x: mcmc_solar.c
	cc mcmc_solar.c -lgsl -lgslcblas -lm -o mcmc_solar.x

clean:
	rm -f *.pdf *.x solar.dat 
