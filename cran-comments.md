
Two bugs were fixed.

Bug #3 in rmnom and rdirmnom led to producing NaN's in output
due to numerical precission issues when the functions were
called with some rather rare combinations of parameters.

Bug #4 was more critical since it led to hanging R when any
of the functions in the package was called with zero-length
input. It was caused by using modulo in implementing
vectorization, where modulo operator in C++ crashes when
calling x % 0. Fix was implemented and all the functions from
now on will work properly and will be tested for this exception.