# Run unit tests
test:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::test()

# CRAN checks
r-cmd-check:
	R CMD check --as-cran .

# Run all examples from the docs
run-examples:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::run_examples()

# Generate test coverage report
coverage:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::test_coverage()

# Run additional checks for the package
check-package:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::check()

# Deploy to CRAN
cran-release: build docs manual _check-for-docs check-package
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::release()

# Build docs
docs:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::document()

_check-for-docs:
	#!/bin/bash
	if [ -z "$(ls -A ./man)" ]; then
		echo "./man directory is empty"
		exit 30
	fi

# Build the PDF manual
manual:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::build_manual()

# Install a local development package
install:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::install()

# Build the package and the manual
build:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::build()
	devtools::build_manual()

# Setup development environment
dev: sys-deps
	#!/usr/bin/env -S Rscript --no-save --no-restore
	install.packages(c('devtools', 'tinytex'))
	devtools::install_deps()
	devtools::install_dev_deps()
	tinytex::install_tinytex()

[linux]
sys-deps:
	sudo apt install \
		libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev libharfbuzz-dev \
		libfribidi-dev libfontconfig1-dev libxml2-dev libssl-dev libcurl4-openssl-dev

# Remove build files
clean:
	rm -rf man/*.Rd
	rm -rf src/*.o
	rm -rf src/*.so
