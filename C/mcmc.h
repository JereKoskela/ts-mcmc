#ifndef MCMC
#define MCMC

#include "tree.h"
#include "recorder.h"
#include <tskit.h>

void kingman_mcmc(recorder *r, tree *t, gsl_rng *gen);

#endif
