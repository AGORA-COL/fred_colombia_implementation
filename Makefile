setup:
	@mkdir -p output
	@mkdir -p scratch
	@echo "Checking FRED..."
	@if [ -d "../FRED" ]; then \
		echo "FRED directory found, skipping clone and make."; \
	else \
		cd .. || (echo "Failed to change directory" && exit 1); \
		git clone https://github.com/confunguido/FRED || (echo "Failed to clone repository" && exit 1); \
		cd FRED/src || (echo "Failed to change into repository directory" && exit 1); \
		make || (echo "Make command failed" && exit 1); \
	fi
	@echo "Checking R version..."
	@R --version | grep 'R version' | cut -d ' ' -f 3 | { read version; if [ "$${version%%.*}" -ge 4 ]; then echo "R version 4 or above is installed"; else echo "R version 4 or above is required"; exit 1; fi; }