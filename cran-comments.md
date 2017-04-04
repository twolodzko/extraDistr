
Two bugs were fixed, including a critical bug that led to
hanging R. Update is provided since the pressing nature of
the second bug.

Bug #3 in rmnom and rdirmnom led to producing NaN's in output
due to numerical precission issues when the functions were
called with some rather rare combinations of parameters.

Bug #4 was more critical since it led to hanging R when any
of the functions in the package was called with zero-length
input. It was caused by using modulo in implementing
vectorization, where modulo operator in C++ crashes when
calling x % 0. Fix was implemented and all the functions from
now on will work properly and will be tested for this exception.

The changes since last update can be compared by following the
link below:
https://github.com/twolodzko/extraDistr/compare/7038ffa41e804a48b0e73a2a85e3e0d69c27f218...master