
solar.pdf: solar.dat 
	python solar.py

solar.dat: mcmc_solar.x
	./solar.x > solar.dat

mcmc_solar.x: mcmc_solar.c
	cc mcmc_solar.c -lm -o mcmc_solar.x

clean:
	rm -f *.pdf *.dat 
