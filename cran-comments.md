## Test environments

* local Ubuntu and Windows, R (release)
* ubuntu 12.04 (on travis-ci), R (oldrel, release, devel)
* Windows (on AppVeyor), R (release)
* win-builder, R (devel)

## R CMD check results

0 errors | 0 warnings | 1 note 
checking installed package size ... NOTE
  installed size is 20.4Mb
  sub-directories of 1Mb or more:
    libs  20.1Mb

## Reverse dependencies

Reverse imports: 	prophet
Reverse suggests: 	greta

The dependencies are not affected. This is bug fix
release and the previous interface of the functions
wa not altered.

## Comments

This is a bug fix release (see NEWS).
