#include "tree.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>

int compare(const void *a, const void *b) {
  int ret = 0;
  if (*(double *)a > *(double *)b) {
    ret = 1;
  } else if (*(double *)a < *(double *)b) {
    ret = -1;
  }
  return ret;
}

int compare_with_index(const void *a, const void *b) {
  int ret = 0;
  sortee *a1 = (sortee *)a;
  sortee *b1 = (sortee *)b;
  if (a1->value > b1->value) {
    ret = 1;
  } else if (a1->value < b1->value) {
    ret = -1;
  }
  return ret;
}

void init_tree(tree *t, const int n) {
  t->left_child = malloc((2 * n - 1) * sizeof(int));
  t->right_child = malloc((2 * n - 1) * sizeof(int));
  t->parent = malloc((2 * n - 1) * sizeof(int));
  t->time = malloc((2 * n - 1) * sizeof(double));
  t->sample_size = n;
  for (int i = 0; i < n; i++) {
    t->left_child[i] = -1;
    t->right_child[i] = -1;
    t->time[i] = 0.0;
  }
  t->parent[2 * n - 2] = -1;
  return;
}

void generate_tree(tree *t, gsl_rng *gen) {
  double r = 0;
  int n = t->sample_size;
  int active_lineages[n];
  for (int i = 0; i < n; i++) {
    active_lineages[i] = i;
  }
  int next_parent = n;
  int l_child, r_child;
  while (n > 1) {
    r += gsl_ran_exponential(gen, 2.0 / (n * (n - 1)));
    l_child = floor(n * gsl_rng_uniform(gen));
    do {
      r_child = floor(n * gsl_rng_uniform(gen));
    } while (l_child == r_child);
    t->parent[active_lineages[l_child]] = next_parent;
    t->parent[active_lineages[r_child]] = next_parent;
    t->left_child[next_parent] = active_lineages[l_child];
    t->right_child[next_parent] = active_lineages[r_child];
    t->time[next_parent] = r;
    active_lineages[(int)fmin(l_child, r_child)] = next_parent;
    for (int i = (int)fmax(l_child, r_child); i < n - 1; i++) {
      active_lineages[i] = active_lineages[i + 1];
    }
    next_parent++;
    n--;
  }
  return;
}

void print_tree(tree *t) {
  for (int i = 0; i < 2 * t->sample_size - 1; i++) {
    printf("%d, %d, %d, %d, %f\n", i, t->parent[i], t->left_child[i],
           t->right_child[i], t->time[i]);
  }
  return;
}

int sibling(tree *t, const int node) {
  int ret = -1;
  int p = t->parent[node];
  if (p != -1) {
    ret = t->left_child[p];
    if (ret == node) {
      ret = t->right_child[p];
    }
  }
  return ret;
}

int grandparent(tree *t, const int node) {
  int ret = -1;
  int p = t->parent[node];
  if (p != -1) {
    ret = t->parent[p];
  }
  return ret;
}

double sample_reattach_time(tree *t, const int child, const int new_sib,
                            gsl_rng *gen) {
  int sib = sibling(t, child);
  int new_sib_parent = t->parent[new_sib];
  int new_sib_grandparent = grandparent(t, new_sib);
  double ret = 0.0;
  if (sib == new_sib) {
    if (new_sib_grandparent == -1) {
      ret = t->time[new_sib] + gsl_ran_exponential(gen, 1.0);
    } else {
      ret =
          t->time[new_sib] + (t->time[new_sib_grandparent] - t->time[new_sib]) *
                                 gsl_rng_uniform(gen);
    }
  } else {
    if (new_sib_parent == -1) {
      ret = t->time[new_sib] + gsl_ran_exponential(gen, 1.0);
    } else {
      ret = t->time[new_sib] +
            (t->time[new_sib_parent] - t->time[new_sib]) * gsl_rng_uniform(gen);
    }
  }
  return ret;
}

double log_reattach_density(tree *t, const int child, const int new_sib,
                            const double new_time) {
  double ret = 0.0;
  int sib = sibling(t, child);
  int new_sib_parent = t->parent[new_sib];
  int new_sib_grandparent = grandparent(t, new_sib);
  if (sib == new_sib) {
    if (new_sib_grandparent == -1) {
      ret = t->time[new_sib] - new_time;
    } else {
      ret = -log(t->time[new_sib_grandparent] - t->time[new_sib]);
    }
  } else {
    if (new_sib_parent == -1) {
      ret = t->time[new_sib] - new_time;
    } else {
      ret = -log(t->time[new_sib_parent] - t->time[new_sib]);
    }
  }
  return ret;
}

double resample_times(tree *t, gsl_rng *gen) {
  double ret = 0.0;
  double lb = 0.0;
  int n = t->sample_size;
  sortee sorted_times[n - 1];
  for (int i = 0; i < n - 1; i++) {
    sorted_times[i].value = t->time[n + i];
    sorted_times[i].index = i;
  }
  qsort(sorted_times, n - 1, sizeof(sortee), compare_with_index);
  double rate = 0.0;
  for (int i = 0; i < n - 1; i++) {
    rate = (n - i) * (n - i - 1) / 2.0;
    t->time[n + sorted_times[i].index] =
        lb + gsl_ran_exponential(gen, 1.0 / rate);
    ret -= rate * (t->time[n + sorted_times[i].index] - lb);
    lb = t->time[n + sorted_times[i].index];
  }
  return ret;
}

double log_resample_times_density(tree *t, const double *r) {
  double ret = 0.0;
  int n = t->sample_size;
  double sorted_times[n];
  for (int i = 0; i < n; i++) {
    sorted_times[i] = r[n + i - 1];
  }
  qsort(sorted_times, n, sizeof(double), compare);
  for (int i = 0; i < n - 1; i++) {
    ret -=
        (sorted_times[i + 1] - sorted_times[i]) * (n - i) * (n - i - 1) / 2.0;
  }
  return ret;
}

int sample_leaf(tree *t, gsl_rng *gen) {
  int n = t->sample_size;
  int ret = floor(n * gsl_rng_uniform(gen));
  return ret;
}

int sample_node(tree *t, gsl_rng *gen) {
  int n = t->sample_size;
  int ret = floor((2 * n - 1) * gsl_rng_uniform(gen));
  return ret;
}

void replace_child(tree *t, const int node, const int child,
                   const int new_child) {
  if (t->left_child[node] == child) {
    t->left_child[node] = new_child;
  } else {
    t->right_child[node] = new_child;
  }
  return;
}

void detach_reattach(tree *t, const int child, const int new_sib,
                     const double new_time) {
  int sib = sibling(t, child);
  int child_parent = t->parent[child];
  int sib_grandparent = grandparent(t, child);
  int new_sib_parent = t->parent[new_sib];
  t->time[child_parent] = new_time;
  if (sib != new_sib) {
    t->parent[new_sib] = child_parent;
    t->parent[sib] = sib_grandparent;
    t->parent[child_parent] = new_sib_parent;
    replace_child(t, child_parent, sib, new_sib);
    if (sib_grandparent != -1) {
      replace_child(t, sib_grandparent, child_parent, sib);
    }
    if (new_sib_parent != -1) {
      replace_child(t, new_sib_parent, new_sib, child_parent);
    }
  }
  return;
}

double log_likelihood(tree *t) {
  double ret = 0.0;
  int n = t->sample_size;
  double sorted_times[n];
  for (int i = 0; i < n; i++) {
    sorted_times[i] = t->time[n + i - 1];
  }
  qsort(sorted_times, n, sizeof(double), compare);
  for (int i = 0; i < n - 1; i++) {
    ret -=
        (n - i) * (n - i - 1) * (sorted_times[i + 1] - sorted_times[i]) / 2.0;
  }
  return ret;
}

void free_tree(tree *t) {
  free(t->left_child);
  free(t->right_child);
  free(t->parent);
  free(t->time);
  free(t);
  return;
}
