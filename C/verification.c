#include "mcmc.h"
#include "recorder.h"
#include "tree.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tskit.h>

#define check_tsk_error(val)                                                   \
  if (val < 0) {                                                               \
    fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));               \
    exit(EXIT_FAILURE);                                                        \
  }

int main() {
  int min_sample = 2;
  int max_sample = 10;
  int steps = 100000;
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));

  double expected_branch_length[max_sample - min_sample + 1];
  expected_branch_length[0] = 2;
  for (int i = min_sample + 1; i <= max_sample; i++) {
    expected_branch_length[i - min_sample] =
        expected_branch_length[i - min_sample - 1] + 2.0 / (i - 1);
  }

  double tmrca, branch_length, clade_prob, tmp, min_time;
  int ret, parent, min_parent;
  for (int i = min_sample; i <= max_sample; i++) {
    tmrca = 0.0;
    branch_length = 0.0;
    clade_prob = 0.0;
    tree *t = malloc(sizeof(tree));
    init_tree(t, i);
    generate_tree(t, gen);
    recorder *r = malloc(sizeof(recorder));
    init_recorder(r, t, steps);
    kingman_mcmc(r, t, gen);
    tree_sequence(r);
    tsk_tree_t tsk_tree;
    ret = tsk_tree_init(&tsk_tree, &r->tree_sequence, 0);
    check_tsk_error(ret);
    for (ret = tsk_tree_first(&tsk_tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tsk_tree)) {
      ret =
          tsk_tree_get_time(&tsk_tree, tsk_tree_get_left_root(&tsk_tree), &tmp);
      check_tsk_error(ret);
      tmrca += tmp * (tsk_tree.interval.right - tsk_tree.interval.left) / steps;
      tsk_tree_get_total_branch_length(&tsk_tree, TSK_NULL, &tmp);
      branch_length +=
          tmp * (tsk_tree.interval.right - tsk_tree.interval.left) / steps;

      ret = tsk_tree_get_parent(&tsk_tree, 0, &parent);
      check_tsk_error(ret);
      ret = tsk_tree_get_time(&tsk_tree, parent, &min_time);
      check_tsk_error(ret);
      min_parent = 0;
      for (int j = 1; j < t->sample_size; j++) {
        ret = tsk_tree_get_parent(&tsk_tree, j, &parent);
        check_tsk_error(ret);
        ret = tsk_tree_get_time(&tsk_tree, parent, &tmp);
        check_tsk_error(ret);
        if (tmp < min_time) {
          min_time = tmp;
          min_parent = j;
        }
      }
      ret = tsk_tree_get_parent(&tsk_tree, min_parent, &min_parent);
      check_tsk_error(ret);
      if ((tsk_tree.left_child[min_parent] == 0 ||
           tsk_tree.left_child[min_parent] == 1) &&
          (tsk_tree.right_child[min_parent] == 0 ||
           tsk_tree.right_child[min_parent] == 1)) {
        clade_prob +=
            (tsk_tree.interval.right - tsk_tree.interval.left) / steps;
      }
    }
    printf("%d %f %f %f %f %f %f\n", i, tmrca, 2 * (1 - 1.0 / i), branch_length,
           expected_branch_length[i - min_sample], log(clade_prob),
           log(2.0 / (i * (i - 1))));
    ret = tsk_tree_free(&tsk_tree);
    check_tsk_error(ret);
    free_recorder(r);
    free_tree(t);
  }
  gsl_rng_free(gen);
  return 0;
}
