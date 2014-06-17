all:      
	cd src && $(MAKE)
	cd src && mv ./bimbam ../bimbam
clean:
	cd src && $(MAKE) clean
	
