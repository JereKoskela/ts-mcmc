#include "mcmc.h"
#include "recorder.h"
#include "tree.h"
#include <tskit.h>

void kingman_mcmc(recorder *r, tree *t, gsl_rng *gen) {
  double llik = log_likelihood(t);
  int child, parent, sib, new_sib;
  double new_time, old_time, alpha;
  double proposal_llik;
  double *old_times = malloc((2 * t->sample_size - 1) * sizeof(double));
  int steps = r->tables.sequence_length;
  for (int i = 1; i < steps; i++) {
    child = sample_leaf(t, gen);
    parent = t->parent[child];
    sib = sibling(t, child);
    new_sib = sib;
    if (t->sample_size > 2) {
      while (new_sib == sib || new_sib == parent || new_sib == child) {
        new_sib = sample_node(t, gen);
      }
    }
    new_time = sample_reattach_time(t, child, new_sib, gen);
    alpha = -log_reattach_density(t, child, new_sib, new_time);
    old_time = t->time[parent];
    detach_reattach(t, child, new_sib, new_time);
    alpha += log_reattach_density(t, child, sib, old_time);
    proposal_llik = log_likelihood(t);
    alpha += proposal_llik - llik;
    if (log(gsl_rng_uniform(gen)) < alpha) {
      llik = proposal_llik;
    } else {
      detach_reattach(t, child, sib, old_time);
    }
    for (int j = t->sample_size; j < 2 * t->sample_size - 1; j++) {
      old_times[j] = t->time[j];
    }
    alpha = -resample_times(t, gen);
    proposal_llik = log_likelihood(t);
    alpha += log_resample_times_density(t, old_times);
    alpha += proposal_llik - llik;
    if (log(gsl_rng_uniform(gen)) < alpha) {
      llik = proposal_llik;
    } else {
      for (int j = t->sample_size; j < 2 * t->sample_size - 1; j++) {
        t->time[j] = old_times[j];
      }
    }
    append_tree(r, t);
  }
  free(old_times);
  return;
}
