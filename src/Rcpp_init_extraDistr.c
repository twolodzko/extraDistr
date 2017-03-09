#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Call calls */
extern SEXP extraDistr_cpp_dbbinom(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbern(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbetapr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbhatt(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbnbinom(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dbpois(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dcat(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddgamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddirichlet(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddirmnom(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddlaplace(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddunif(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ddweibull(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dfatigue(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dfrechet(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dgev(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dgompertz(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dgpd(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dgpois(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dgumbel(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dhcauchy(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dhnorm(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dht(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dhuber(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dinvgamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dkumar(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dlaplace(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dlgser(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dlomax(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dmixnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dmixpois(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dmnom(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dmvhyper(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dnhyper(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dnsbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dnst(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dpareto(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dpower(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dprop(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_drayleigh(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dsgomp(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dskellam(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dslash(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dtbinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dtnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dtpois(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dtriang(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dwald(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dzib(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dzinb(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_dzip(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pbbinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pbern(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pbetapr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pbhatt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pbnbinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pcat(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pdlaplace(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pdunif(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pdweibull(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pfatigue(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pfrechet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pgev(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pgompertz(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pgpd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pgpois(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pgumbel(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_phcauchy(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_phnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pht(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_phuber(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pkumar(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_plaplace(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_plgser(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_plomax(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pmixnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pmixpois(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pnhyper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pnsbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pnst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ppareto(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ppower(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pprop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_prayleigh(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_psgomp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pslash(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ptbinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ptnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ptpois(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_ptriang(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pwald(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pzib(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pzinb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_pzip(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qbern(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qbetapr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qcat(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qdunif(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qdweibull(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qfatigue(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qfrechet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qgev(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qgompertz(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qgpd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qgumbel(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qhcauchy(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qhnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qht(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qhuber(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qkumar(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qlaplace(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qlgser(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qlomax(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qnhyper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qnsbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qnst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qpareto(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qpower(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qprop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qrayleigh(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qtbinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qtlambda(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qtnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qtpois(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qtriang(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qzib(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qzinb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_qzip(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbern(SEXP, SEXP);
extern SEXP extraDistr_cpp_rbetapr(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbhatt(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbnbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rbpois(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rcat(SEXP, SEXP);
extern SEXP extraDistr_cpp_rcatlp(SEXP, SEXP);
extern SEXP extraDistr_cpp_rdirichlet(SEXP, SEXP);
extern SEXP extraDistr_cpp_rdirmnom(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rdlaplace(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rdunif(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rdweibull(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rfatigue(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rfrechet(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rgev(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rgompertz(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rgpd(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rgpois(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rgumbel(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rhcauchy(SEXP, SEXP);
extern SEXP extraDistr_cpp_rhnorm(SEXP, SEXP);
extern SEXP extraDistr_cpp_rht(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rhuber(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rkumar(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rlaplace(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rlgser(SEXP, SEXP);
extern SEXP extraDistr_cpp_rlomax(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rmixnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rmixpois(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rmnom(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rmvhyper(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rnhyper(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rnsbeta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rnst(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rpareto(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rpower(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rprop(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rrayleigh(SEXP, SEXP);
extern SEXP extraDistr_cpp_rsgomp(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rsign(SEXP);
extern SEXP extraDistr_cpp_rskellam(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rslash(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rtbinom(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rtlambda(SEXP, SEXP);
extern SEXP extraDistr_cpp_rtnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rtpois(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rtriang(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rwald(SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rzib(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rzinb(SEXP, SEXP, SEXP, SEXP);
extern SEXP extraDistr_cpp_rzip(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"extraDistr_cpp_dbbinom",    (DL_FUNC) &extraDistr_cpp_dbbinom,    5},
  {"extraDistr_cpp_dbern",      (DL_FUNC) &extraDistr_cpp_dbern,      3},
  {"extraDistr_cpp_dbetapr",    (DL_FUNC) &extraDistr_cpp_dbetapr,    5},
  {"extraDistr_cpp_dbhatt",     (DL_FUNC) &extraDistr_cpp_dbhatt,     5},
  {"extraDistr_cpp_dbnbinom",   (DL_FUNC) &extraDistr_cpp_dbnbinom,   5},
  {"extraDistr_cpp_dbnorm",     (DL_FUNC) &extraDistr_cpp_dbnorm,     8},
  {"extraDistr_cpp_dbpois",     (DL_FUNC) &extraDistr_cpp_dbpois,     6},
  {"extraDistr_cpp_dcat",       (DL_FUNC) &extraDistr_cpp_dcat,       3},
  {"extraDistr_cpp_ddgamma",    (DL_FUNC) &extraDistr_cpp_ddgamma,    4},
  {"extraDistr_cpp_ddirichlet", (DL_FUNC) &extraDistr_cpp_ddirichlet, 3},
  {"extraDistr_cpp_ddirmnom",   (DL_FUNC) &extraDistr_cpp_ddirmnom,   4},
  {"extraDistr_cpp_ddlaplace",  (DL_FUNC) &extraDistr_cpp_ddlaplace,  4},
  {"extraDistr_cpp_ddnorm",     (DL_FUNC) &extraDistr_cpp_ddnorm,     4},
  {"extraDistr_cpp_ddunif",     (DL_FUNC) &extraDistr_cpp_ddunif,     4},
  {"extraDistr_cpp_ddweibull",  (DL_FUNC) &extraDistr_cpp_ddweibull,  4},
  {"extraDistr_cpp_dfatigue",   (DL_FUNC) &extraDistr_cpp_dfatigue,   5},
  {"extraDistr_cpp_dfrechet",   (DL_FUNC) &extraDistr_cpp_dfrechet,   5},
  {"extraDistr_cpp_dgev",       (DL_FUNC) &extraDistr_cpp_dgev,       5},
  {"extraDistr_cpp_dgompertz",  (DL_FUNC) &extraDistr_cpp_dgompertz,  4},
  {"extraDistr_cpp_dgpd",       (DL_FUNC) &extraDistr_cpp_dgpd,       5},
  {"extraDistr_cpp_dgpois",     (DL_FUNC) &extraDistr_cpp_dgpois,     4},
  {"extraDistr_cpp_dgumbel",    (DL_FUNC) &extraDistr_cpp_dgumbel,    4},
  {"extraDistr_cpp_dhcauchy",   (DL_FUNC) &extraDistr_cpp_dhcauchy,   3},
  {"extraDistr_cpp_dhnorm",     (DL_FUNC) &extraDistr_cpp_dhnorm,     3},
  {"extraDistr_cpp_dht",        (DL_FUNC) &extraDistr_cpp_dht,        4},
  {"extraDistr_cpp_dhuber",     (DL_FUNC) &extraDistr_cpp_dhuber,     5},
  {"extraDistr_cpp_dinvgamma",  (DL_FUNC) &extraDistr_cpp_dinvgamma,  4},
  {"extraDistr_cpp_dkumar",     (DL_FUNC) &extraDistr_cpp_dkumar,     4},
  {"extraDistr_cpp_dlaplace",   (DL_FUNC) &extraDistr_cpp_dlaplace,   4},
  {"extraDistr_cpp_dlgser",     (DL_FUNC) &extraDistr_cpp_dlgser,     3},
  {"extraDistr_cpp_dlomax",     (DL_FUNC) &extraDistr_cpp_dlomax,     4},
  {"extraDistr_cpp_dmixnorm",   (DL_FUNC) &extraDistr_cpp_dmixnorm,   5},
  {"extraDistr_cpp_dmixpois",   (DL_FUNC) &extraDistr_cpp_dmixpois,   4},
  {"extraDistr_cpp_dmnom",      (DL_FUNC) &extraDistr_cpp_dmnom,      4},
  {"extraDistr_cpp_dmvhyper",   (DL_FUNC) &extraDistr_cpp_dmvhyper,   4},
  {"extraDistr_cpp_dnhyper",    (DL_FUNC) &extraDistr_cpp_dnhyper,    5},
  {"extraDistr_cpp_dnsbeta",    (DL_FUNC) &extraDistr_cpp_dnsbeta,    6},
  {"extraDistr_cpp_dnst",       (DL_FUNC) &extraDistr_cpp_dnst,       5},
  {"extraDistr_cpp_dpareto",    (DL_FUNC) &extraDistr_cpp_dpareto,    4},
  {"extraDistr_cpp_dpower",     (DL_FUNC) &extraDistr_cpp_dpower,     4},
  {"extraDistr_cpp_dprop",      (DL_FUNC) &extraDistr_cpp_dprop,      5},
  {"extraDistr_cpp_drayleigh",  (DL_FUNC) &extraDistr_cpp_drayleigh,  3},
  {"extraDistr_cpp_dsgomp",     (DL_FUNC) &extraDistr_cpp_dsgomp,     4},
  {"extraDistr_cpp_dskellam",   (DL_FUNC) &extraDistr_cpp_dskellam,   4},
  {"extraDistr_cpp_dslash",     (DL_FUNC) &extraDistr_cpp_dslash,     4},
  {"extraDistr_cpp_dtbinom",    (DL_FUNC) &extraDistr_cpp_dtbinom,    6},
  {"extraDistr_cpp_dtnorm",     (DL_FUNC) &extraDistr_cpp_dtnorm,     6},
  {"extraDistr_cpp_dtpois",     (DL_FUNC) &extraDistr_cpp_dtpois,     5},
  {"extraDistr_cpp_dtriang",    (DL_FUNC) &extraDistr_cpp_dtriang,    5},
  {"extraDistr_cpp_dwald",      (DL_FUNC) &extraDistr_cpp_dwald,      4},
  {"extraDistr_cpp_dzib",       (DL_FUNC) &extraDistr_cpp_dzib,       5},
  {"extraDistr_cpp_dzinb",      (DL_FUNC) &extraDistr_cpp_dzinb,      5},
  {"extraDistr_cpp_dzip",       (DL_FUNC) &extraDistr_cpp_dzip,       4},
  {"extraDistr_cpp_pbbinom",    (DL_FUNC) &extraDistr_cpp_pbbinom,    6},
  {"extraDistr_cpp_pbern",      (DL_FUNC) &extraDistr_cpp_pbern,      4},
  {"extraDistr_cpp_pbetapr",    (DL_FUNC) &extraDistr_cpp_pbetapr,    6},
  {"extraDistr_cpp_pbhatt",     (DL_FUNC) &extraDistr_cpp_pbhatt,     6},
  {"extraDistr_cpp_pbnbinom",   (DL_FUNC) &extraDistr_cpp_pbnbinom,   6},
  {"extraDistr_cpp_pcat",       (DL_FUNC) &extraDistr_cpp_pcat,       4},
  {"extraDistr_cpp_pdlaplace",  (DL_FUNC) &extraDistr_cpp_pdlaplace,  5},
  {"extraDistr_cpp_pdunif",     (DL_FUNC) &extraDistr_cpp_pdunif,     5},
  {"extraDistr_cpp_pdweibull",  (DL_FUNC) &extraDistr_cpp_pdweibull,  5},
  {"extraDistr_cpp_pfatigue",   (DL_FUNC) &extraDistr_cpp_pfatigue,   6},
  {"extraDistr_cpp_pfrechet",   (DL_FUNC) &extraDistr_cpp_pfrechet,   6},
  {"extraDistr_cpp_pgev",       (DL_FUNC) &extraDistr_cpp_pgev,       6},
  {"extraDistr_cpp_pgompertz",  (DL_FUNC) &extraDistr_cpp_pgompertz,  5},
  {"extraDistr_cpp_pgpd",       (DL_FUNC) &extraDistr_cpp_pgpd,       6},
  {"extraDistr_cpp_pgpois",     (DL_FUNC) &extraDistr_cpp_pgpois,     5},
  {"extraDistr_cpp_pgumbel",    (DL_FUNC) &extraDistr_cpp_pgumbel,    5},
  {"extraDistr_cpp_phcauchy",   (DL_FUNC) &extraDistr_cpp_phcauchy,   4},
  {"extraDistr_cpp_phnorm",     (DL_FUNC) &extraDistr_cpp_phnorm,     4},
  {"extraDistr_cpp_pht",        (DL_FUNC) &extraDistr_cpp_pht,        5},
  {"extraDistr_cpp_phuber",     (DL_FUNC) &extraDistr_cpp_phuber,     6},
  {"extraDistr_cpp_pkumar",     (DL_FUNC) &extraDistr_cpp_pkumar,     5},
  {"extraDistr_cpp_plaplace",   (DL_FUNC) &extraDistr_cpp_plaplace,   5},
  {"extraDistr_cpp_plgser",     (DL_FUNC) &extraDistr_cpp_plgser,     4},
  {"extraDistr_cpp_plomax",     (DL_FUNC) &extraDistr_cpp_plomax,     5},
  {"extraDistr_cpp_pmixnorm",   (DL_FUNC) &extraDistr_cpp_pmixnorm,   6},
  {"extraDistr_cpp_pmixpois",   (DL_FUNC) &extraDistr_cpp_pmixpois,   5},
  {"extraDistr_cpp_pnhyper",    (DL_FUNC) &extraDistr_cpp_pnhyper,    6},
  {"extraDistr_cpp_pnsbeta",    (DL_FUNC) &extraDistr_cpp_pnsbeta,    7},
  {"extraDistr_cpp_pnst",       (DL_FUNC) &extraDistr_cpp_pnst,       6},
  {"extraDistr_cpp_ppareto",    (DL_FUNC) &extraDistr_cpp_ppareto,    5},
  {"extraDistr_cpp_ppower",     (DL_FUNC) &extraDistr_cpp_ppower,     5},
  {"extraDistr_cpp_pprop",      (DL_FUNC) &extraDistr_cpp_pprop,      6},
  {"extraDistr_cpp_prayleigh",  (DL_FUNC) &extraDistr_cpp_prayleigh,  4},
  {"extraDistr_cpp_psgomp",     (DL_FUNC) &extraDistr_cpp_psgomp,     5},
  {"extraDistr_cpp_pslash",     (DL_FUNC) &extraDistr_cpp_pslash,     5},
  {"extraDistr_cpp_ptbinom",    (DL_FUNC) &extraDistr_cpp_ptbinom,    7},
  {"extraDistr_cpp_ptnorm",     (DL_FUNC) &extraDistr_cpp_ptnorm,     7},
  {"extraDistr_cpp_ptpois",     (DL_FUNC) &extraDistr_cpp_ptpois,     6},
  {"extraDistr_cpp_ptriang",    (DL_FUNC) &extraDistr_cpp_ptriang,    6},
  {"extraDistr_cpp_pwald",      (DL_FUNC) &extraDistr_cpp_pwald,      5},
  {"extraDistr_cpp_pzib",       (DL_FUNC) &extraDistr_cpp_pzib,       6},
  {"extraDistr_cpp_pzinb",      (DL_FUNC) &extraDistr_cpp_pzinb,      6},
  {"extraDistr_cpp_pzip",       (DL_FUNC) &extraDistr_cpp_pzip,       5},
  {"extraDistr_cpp_qbern",      (DL_FUNC) &extraDistr_cpp_qbern,      4},
  {"extraDistr_cpp_qbetapr",    (DL_FUNC) &extraDistr_cpp_qbetapr,    6},
  {"extraDistr_cpp_qcat",       (DL_FUNC) &extraDistr_cpp_qcat,       4},
  {"extraDistr_cpp_qdunif",     (DL_FUNC) &extraDistr_cpp_qdunif,     5},
  {"extraDistr_cpp_qdweibull",  (DL_FUNC) &extraDistr_cpp_qdweibull,  5},
  {"extraDistr_cpp_qfatigue",   (DL_FUNC) &extraDistr_cpp_qfatigue,   6},
  {"extraDistr_cpp_qfrechet",   (DL_FUNC) &extraDistr_cpp_qfrechet,   6},
  {"extraDistr_cpp_qgev",       (DL_FUNC) &extraDistr_cpp_qgev,       6},
  {"extraDistr_cpp_qgompertz",  (DL_FUNC) &extraDistr_cpp_qgompertz,  5},
  {"extraDistr_cpp_qgpd",       (DL_FUNC) &extraDistr_cpp_qgpd,       6},
  {"extraDistr_cpp_qgumbel",    (DL_FUNC) &extraDistr_cpp_qgumbel,    5},
  {"extraDistr_cpp_qhcauchy",   (DL_FUNC) &extraDistr_cpp_qhcauchy,   4},
  {"extraDistr_cpp_qhnorm",     (DL_FUNC) &extraDistr_cpp_qhnorm,     4},
  {"extraDistr_cpp_qht",        (DL_FUNC) &extraDistr_cpp_qht,        5},
  {"extraDistr_cpp_qhuber",     (DL_FUNC) &extraDistr_cpp_qhuber,     6},
  {"extraDistr_cpp_qkumar",     (DL_FUNC) &extraDistr_cpp_qkumar,     5},
  {"extraDistr_cpp_qlaplace",   (DL_FUNC) &extraDistr_cpp_qlaplace,   5},
  {"extraDistr_cpp_qlgser",     (DL_FUNC) &extraDistr_cpp_qlgser,     4},
  {"extraDistr_cpp_qlomax",     (DL_FUNC) &extraDistr_cpp_qlomax,     5},
  {"extraDistr_cpp_qnhyper",    (DL_FUNC) &extraDistr_cpp_qnhyper,    6},
  {"extraDistr_cpp_qnsbeta",    (DL_FUNC) &extraDistr_cpp_qnsbeta,    7},
  {"extraDistr_cpp_qnst",       (DL_FUNC) &extraDistr_cpp_qnst,       6},
  {"extraDistr_cpp_qpareto",    (DL_FUNC) &extraDistr_cpp_qpareto,    5},
  {"extraDistr_cpp_qpower",     (DL_FUNC) &extraDistr_cpp_qpower,     5},
  {"extraDistr_cpp_qprop",      (DL_FUNC) &extraDistr_cpp_qprop,      6},
  {"extraDistr_cpp_qrayleigh",  (DL_FUNC) &extraDistr_cpp_qrayleigh,  4},
  {"extraDistr_cpp_qtbinom",    (DL_FUNC) &extraDistr_cpp_qtbinom,    7},
  {"extraDistr_cpp_qtlambda",   (DL_FUNC) &extraDistr_cpp_qtlambda,   4},
  {"extraDistr_cpp_qtnorm",     (DL_FUNC) &extraDistr_cpp_qtnorm,     7},
  {"extraDistr_cpp_qtpois",     (DL_FUNC) &extraDistr_cpp_qtpois,     6},
  {"extraDistr_cpp_qtriang",    (DL_FUNC) &extraDistr_cpp_qtriang,    6},
  {"extraDistr_cpp_qzib",       (DL_FUNC) &extraDistr_cpp_qzib,       6},
  {"extraDistr_cpp_qzinb",      (DL_FUNC) &extraDistr_cpp_qzinb,      6},
  {"extraDistr_cpp_qzip",       (DL_FUNC) &extraDistr_cpp_qzip,       5},
  {"extraDistr_cpp_rbbinom",    (DL_FUNC) &extraDistr_cpp_rbbinom,    4},
  {"extraDistr_cpp_rbern",      (DL_FUNC) &extraDistr_cpp_rbern,      2},
  {"extraDistr_cpp_rbetapr",    (DL_FUNC) &extraDistr_cpp_rbetapr,    4},
  {"extraDistr_cpp_rbhatt",     (DL_FUNC) &extraDistr_cpp_rbhatt,     4},
  {"extraDistr_cpp_rbnbinom",   (DL_FUNC) &extraDistr_cpp_rbnbinom,   4},
  {"extraDistr_cpp_rbnorm",     (DL_FUNC) &extraDistr_cpp_rbnorm,     6},
  {"extraDistr_cpp_rbpois",     (DL_FUNC) &extraDistr_cpp_rbpois,     4},
  {"extraDistr_cpp_rcat",       (DL_FUNC) &extraDistr_cpp_rcat,       2},
  {"extraDistr_cpp_rcatlp",     (DL_FUNC) &extraDistr_cpp_rcatlp,     2},
  {"extraDistr_cpp_rdirichlet", (DL_FUNC) &extraDistr_cpp_rdirichlet, 2},
  {"extraDistr_cpp_rdirmnom",   (DL_FUNC) &extraDistr_cpp_rdirmnom,   3},
  {"extraDistr_cpp_rdlaplace",  (DL_FUNC) &extraDistr_cpp_rdlaplace,  3},
  {"extraDistr_cpp_rdunif",     (DL_FUNC) &extraDistr_cpp_rdunif,     3},
  {"extraDistr_cpp_rdweibull",  (DL_FUNC) &extraDistr_cpp_rdweibull,  3},
  {"extraDistr_cpp_rfatigue",   (DL_FUNC) &extraDistr_cpp_rfatigue,   4},
  {"extraDistr_cpp_rfrechet",   (DL_FUNC) &extraDistr_cpp_rfrechet,   4},
  {"extraDistr_cpp_rgev",       (DL_FUNC) &extraDistr_cpp_rgev,       4},
  {"extraDistr_cpp_rgompertz",  (DL_FUNC) &extraDistr_cpp_rgompertz,  3},
  {"extraDistr_cpp_rgpd",       (DL_FUNC) &extraDistr_cpp_rgpd,       4},
  {"extraDistr_cpp_rgpois",     (DL_FUNC) &extraDistr_cpp_rgpois,     3},
  {"extraDistr_cpp_rgumbel",    (DL_FUNC) &extraDistr_cpp_rgumbel,    3},
  {"extraDistr_cpp_rhcauchy",   (DL_FUNC) &extraDistr_cpp_rhcauchy,   2},
  {"extraDistr_cpp_rhnorm",     (DL_FUNC) &extraDistr_cpp_rhnorm,     2},
  {"extraDistr_cpp_rht",        (DL_FUNC) &extraDistr_cpp_rht,        3},
  {"extraDistr_cpp_rhuber",     (DL_FUNC) &extraDistr_cpp_rhuber,     4},
  {"extraDistr_cpp_rkumar",     (DL_FUNC) &extraDistr_cpp_rkumar,     3},
  {"extraDistr_cpp_rlaplace",   (DL_FUNC) &extraDistr_cpp_rlaplace,   3},
  {"extraDistr_cpp_rlgser",     (DL_FUNC) &extraDistr_cpp_rlgser,     2},
  {"extraDistr_cpp_rlomax",     (DL_FUNC) &extraDistr_cpp_rlomax,     3},
  {"extraDistr_cpp_rmixnorm",   (DL_FUNC) &extraDistr_cpp_rmixnorm,   4},
  {"extraDistr_cpp_rmixpois",   (DL_FUNC) &extraDistr_cpp_rmixpois,   3},
  {"extraDistr_cpp_rmnom",      (DL_FUNC) &extraDistr_cpp_rmnom,      3},
  {"extraDistr_cpp_rmvhyper",   (DL_FUNC) &extraDistr_cpp_rmvhyper,   3},
  {"extraDistr_cpp_rnhyper",    (DL_FUNC) &extraDistr_cpp_rnhyper,    4},
  {"extraDistr_cpp_rnsbeta",    (DL_FUNC) &extraDistr_cpp_rnsbeta,    5},
  {"extraDistr_cpp_rnst",       (DL_FUNC) &extraDistr_cpp_rnst,       4},
  {"extraDistr_cpp_rpareto",    (DL_FUNC) &extraDistr_cpp_rpareto,    3},
  {"extraDistr_cpp_rpower",     (DL_FUNC) &extraDistr_cpp_rpower,     3},
  {"extraDistr_cpp_rprop",      (DL_FUNC) &extraDistr_cpp_rprop,      4},
  {"extraDistr_cpp_rrayleigh",  (DL_FUNC) &extraDistr_cpp_rrayleigh,  2},
  {"extraDistr_cpp_rsgomp",     (DL_FUNC) &extraDistr_cpp_rsgomp,     3},
  {"extraDistr_cpp_rsign",      (DL_FUNC) &extraDistr_cpp_rsign,      1},
  {"extraDistr_cpp_rskellam",   (DL_FUNC) &extraDistr_cpp_rskellam,   3},
  {"extraDistr_cpp_rslash",     (DL_FUNC) &extraDistr_cpp_rslash,     3},
  {"extraDistr_cpp_rtbinom",    (DL_FUNC) &extraDistr_cpp_rtbinom,    5},
  {"extraDistr_cpp_rtlambda",   (DL_FUNC) &extraDistr_cpp_rtlambda,   2},
  {"extraDistr_cpp_rtnorm",     (DL_FUNC) &extraDistr_cpp_rtnorm,     5},
  {"extraDistr_cpp_rtpois",     (DL_FUNC) &extraDistr_cpp_rtpois,     4},
  {"extraDistr_cpp_rtriang",    (DL_FUNC) &extraDistr_cpp_rtriang,    4},
  {"extraDistr_cpp_rwald",      (DL_FUNC) &extraDistr_cpp_rwald,      3},
  {"extraDistr_cpp_rzib",       (DL_FUNC) &extraDistr_cpp_rzib,       4},
  {"extraDistr_cpp_rzinb",      (DL_FUNC) &extraDistr_cpp_rzinb,      4},
  {"extraDistr_cpp_rzip",       (DL_FUNC) &extraDistr_cpp_rzip,       3},
  {NULL, NULL, 0}
};

void R_init_extraDistr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}