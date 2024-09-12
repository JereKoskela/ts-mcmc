#ifndef TREE
#define TREE

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

typedef struct {
  double value;
  int index;
} sortee;

typedef struct {
  int *left_child;
  int *right_child;
  int *parent;
  double *time;
  int sample_size;
} tree;

int compare(const void *a, const void *b);
int compare_with_index(const void *a, const void *b);

void init_tree(tree *t, const int n);
void generate_tree(tree *t, gsl_rng *gen);
void print_tree(tree *t);
int sibling(tree *t, const int node);
int grandparent(tree *t, const int node);
double sample_reattach_time(tree *t, const int child, const int new_sib,
                            gsl_rng *gen);
double log_reattach_density(tree *t, const int child, const int new_sib,
                            const double new_time);
double resample_times(tree *t, gsl_rng *gen);
double log_resample_times_density(tree *t, const double *r);
int sample_leaf(tree *t, gsl_rng *gen);
int sample_node(tree *t, gsl_rng *gen);
void replace_child(tree *t, const int node, const int child,
                   const int new_child);
void detach_reattach(tree *t, const int child, const int new_sib,
                     const double new_time);
double log_likelihood(tree *t);
void free_tree(tree *t);

#endif
