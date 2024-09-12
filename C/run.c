#include "mcmc.h"
#include "recorder.h"
#include "tree.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tskit.h>

int main() {
  int sample_size = 10;
  int steps = 10;

  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));

  tree *t = malloc(sizeof(tree));
  init_tree(t, sample_size);
  generate_tree(t, gen);

  recorder *r = malloc(sizeof(recorder));
  init_recorder(r, t, steps);

  kingman_mcmc(r, t, gen);
  tree_sequence(r);
  tsk_treeseq_dump(&r->tree_sequence, "my_ts.trees", 0);

  free_tree(t);
  free_recorder(r);
  gsl_rng_free(gen);

  return 0;
}
