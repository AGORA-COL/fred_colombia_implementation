setup:
	@mkdir -p output
	@mkdir -p scratch
	@mkdir -p fred_run_stages/run_files
	@echo "Checking FRED..."
	@if [ -d "../FRED" ]; then \
		echo "FRED directory found, skipping clone and make. Remember to add the synthetic populations"; \
	else \
		cd .. || (echo "Failed to change directory" && exit 1); \
		git clone https://github.com/confunguido/FRED || (echo "Failed to clone repository" && exit 1); \
        echo "Copying custom Epidemic.cc and Global.h files"; \
        cp -R ./fred_colombia_implementation/input_files/FRED_compile_files/* FRED/src/ || (echo "Failed to copy files" && exit 1); \
		cd FRED/src || (echo "Failed to change into repository directory" && exit 1); \
		make || (echo "Make command failed" && exit 1); \
		echo "Remember to add the synthetic populations"; \
	fi
	@echo "Checking R version..."
	@R --version | grep 'R version' | cut -d ' ' -f 3 | { read version; if [ "$${version%%.*}" -ge 4 ]; then echo "R version 4 or above is installed"; else echo "R version 4 or above is required"; exit 1; fi; }

PATH_TO_EXTRACT := ../FRED/populations/
download:
	@output=$$(python scripts/synthetic_populations_setup.py "$(files)" "$(PATH_TO_EXTRACT)"); \
	downloaded_files=$$(echo $$output | cut -d'|' -f1); \
	skipped_files=$$(echo $$output | cut -d'|' -f2); \
	for file in $$downloaded_files; do \
		echo "Extracting $$file"; \
		unzip -o $$file -d $(PATH_TO_EXTRACT)`basename $$file .zip`; \
		echo "Deleting $$file.zip"; \
		rm $$file.zip; \
	done; \
	if [ "$$skipped_files" ]; then \
		echo "Skipped files: $$skipped_files"; \
	fi