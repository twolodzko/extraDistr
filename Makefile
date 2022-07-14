
.PHONY: test
test: ## Run unit tests
	R --vanilla --silent -e "devtools::test()"

.PHONY: coverage
coverage: ## Generate test coverage report
	R --vanilla --silent -e "devtools::test_coverage()"

.PHONY: check-pacakge
check-pacakge: ## Run additional checks for the package
	R --vanilla --silent -e "devtools::check()"

.PHONY: cran-release
cran-release: ## Deploy to CRAN
	R --vanilla --silent -e "devtools::release()"

.PHONY: install
install: ## Install a local development package
	R --silent --slave --no-save --no-restore -e "devtools::install()"

.PHONY: build
build: ## Build the package and the manual
	R --vanilla -e "devtools::build()"
	R --vanilla -e "devtools::build_manual()"

.PHONY: dev
dev: ## Setup development environment
	R --silent --slave --no-save --no-restore -e "install.packages('devtools')"
	R --silent --slave --no-save --no-restore -e "devtools::install_deps()"
	R --silent --slave --no-save --no-restore -e "devtools::install_dev_deps() "

.PHONY: help
help: ## Print help
	@ grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
