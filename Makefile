setup:
	mkdir output_files
	mkdir scratch
	cd .. && \
	git clone https://github.com/confunguido/FRED
	cd FRED/src && \
	make
