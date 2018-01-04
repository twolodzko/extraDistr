## Test environments

* local Ubuntu and Windows, R (release)
* ubuntu 12.04 (on travis-ci), R (oldrel, release, devel)
* Windows (on AppVeyor), R (release)
* win-builder, R (devel)

## R CMD check results

0 errors | 0 warnings | 1 note 
checking installed package size ... NOTE
  installed size is 20.7Mb
  sub-directories of 1Mb or more:
    libs  20.3Mb

## Reverse dependencies

Reverse depends: 	vfcp
Reverse imports: 	prophet, trend
Reverse suggests: 	greta, SimDesign

The dependencies should not be affected.

The revdep tests on greta resulted in warnings and errors,
but to the best of my knowledge, this is not due to changes
in my package. Greta uses my package as a benchmark
for testing implementation of few of the distributions.
In the current version of my package I added more tests, 
including tests comparing my package to naive pure-R
implementations and tested it by comparing it to 
other popular packages (e.g. VGAM, evd). The functions
from greta package that were tested against functions
in my package were tested against the VGAM implementations
and didn't lead to any errors.

## Comments

This is a bug fix release (see NEWS).

The deprecated functions from pervious version of the package
were removed in this version.


