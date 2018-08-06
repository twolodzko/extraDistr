## Test environments

* local Ubuntu and Windows, R (release)
* ubuntu 12.04 (on travis-ci), R (oldrel, release, devel)
* Windows (on AppVeyor), R (release)
* win-builder, R (devel)

## R CMD check results

0 errors | 0 warnings | 1 note 
checking installed package size ... NOTE
  installed size is 21.0Mb
  sub-directories of 1Mb or more:
    libs  20.6Mb

## Reverse dependencies

Reverse depends: 	vfcp
Reverse imports: 	prophet, trend
Reverse suggests: 	greta, SimDesign

The dependencies should not be affected.

## Comments

Updated the DESCRPTION file to mention packages used in
unit tests in "Suggests". Fixed a bug in one function.


