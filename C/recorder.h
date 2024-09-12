#ifndef RECORDER
#define RECORDER

#include "tree.h"
#include <tskit.h>

typedef struct {
  tsk_table_collection_t tables;
  tsk_treeseq_t tree_sequence;
  int num_nodes;
  int *tree_to_table_map;
  int *left_edge_below;
  int *right_edge_below;
} recorder;

void init_recorder(recorder *r, tree *t, const int len);
void append_tree(recorder *r, tree *t);
void tree_sequence(recorder *r);
void free_recorder(recorder *r);

#endif
