# Run unit tests
test:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::test()

# Run all examples from the docs
run-examples:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::run_examples()

# Generate test coverage report
coverage:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::test_coverage()

# Run additional checks for the package
check-pacakge:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::check()

# Deploy to CRAN
cran-release:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::release()

# Build docs
document:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	devtools::document()
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
dev:
	#!/usr/bin/env -S Rscript --no-save --no-restore
	install.packages('devtools')
	devtools::install_deps()
	devtools::install_dev_deps()

clean:
	rm -rf man/*.Rd
	rm -rf src/*.o
	rm -rf src/*.so
