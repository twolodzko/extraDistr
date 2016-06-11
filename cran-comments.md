
## Update

Namespace was cleaned, all math functions are called from
std library. Rccp namespace is not used globally any more.
All numeric values are now explicitely <double> or <int>
(or casted to) to avoid overloading ambiguity for math functions.

Moreover, C++ code was cleaned up. There were some improvements
in several functions, e.g. for mixture of normals, or mixture of
Poissons, discrete uniform, or categorical distribution.

Vignette with full list of implemented distributions was added.