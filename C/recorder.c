#include "recorder.h"
#include "tree.h"
#include <tskit.h>

#define check_tsk_error(val)                                                   \
  if (val < 0) {                                                               \
    fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));               \
    exit(EXIT_FAILURE);                                                        \
  }

void init_recorder(recorder *r, tree *t, const int len) {
  int ret;
  ret = tsk_table_collection_init(&r->tables, 0);
  check_tsk_error(ret);
  r->tables.sequence_length = len;
  r->num_nodes = 2 * t->sample_size - 1;
  r->tree_to_table_map = malloc(r->num_nodes * sizeof(int));
  r->left_edge_below = malloc(r->num_nodes * sizeof(int));
  r->right_edge_below = malloc(r->num_nodes * sizeof(int));
  for (int i = 0; i < t->sample_size; i++) {
    ret = tsk_node_table_add_row(&r->tables.nodes, TSK_NODE_IS_SAMPLE, 0,
                                 TSK_NULL, TSK_NULL, NULL, 0);
    check_tsk_error(ret);
    r->tree_to_table_map[i] = i;
    r->left_edge_below[i] = -1;
    r->right_edge_below[i] = -1;
  }
  for (int i = t->sample_size; i < r->num_nodes; i++) {
    ret = tsk_node_table_add_row(&r->tables.nodes, 0, t->time[i], TSK_NULL,
                                 TSK_NULL, NULL, 0);
    check_tsk_error(ret);
    r->tree_to_table_map[i] = i;
    ret = tsk_edge_table_add_row(&r->tables.edges, 0, 1, i, t->left_child[i],
                                 NULL, 0);
    check_tsk_error(ret);
    r->left_edge_below[i] = r->tables.edges.num_rows - 1;
    ret = tsk_edge_table_add_row(&r->tables.edges, 0, 1, i, t->right_child[i],
                                 NULL, 0);
    check_tsk_error(ret);
    r->right_edge_below[i] = r->tables.edges.num_rows - 1;
  }
  return;
}

void append_tree(recorder *r, tree *t) {
  int ret;
  int right = r->tables.edges.right[r->left_edge_below[r->num_nodes - 1]];
  for (int i = t->sample_size; i < 2 * t->sample_size - 1; i++) {
    // We assume MCMC tree updates always change node time
    // False positives from floating point comparison inflate memory cost a bit
    if (t->time[i] != r->tables.nodes.time[r->tree_to_table_map[i]]) {
      ret = tsk_node_table_add_row(&r->tables.nodes, 0, t->time[i], TSK_NULL,
                                   TSK_NULL, NULL, 0);
      check_tsk_error(ret);
      r->tree_to_table_map[i] = r->tables.nodes.num_rows - 1;
    }
  }
  int tree_parent, tree_child;
  int edge, edge_child, edge_parent;
  for (int i = t->sample_size; i < 2 * t->sample_size - 1; i++) {
    tree_parent = r->tree_to_table_map[i];
    tree_child = r->tree_to_table_map[t->left_child[i]];
    edge = r->left_edge_below[i];
    edge_child = r->tables.edges.child[edge];
    edge_parent = r->tables.edges.parent[edge];
    if (tree_child != edge_child || tree_parent != edge_parent) {
      ret = tsk_edge_table_add_row(&r->tables.edges, right, right + 1,
                                   tree_parent, tree_child, NULL, 0);
      check_tsk_error(ret);
      r->left_edge_below[i] = r->tables.edges.num_rows - 1;
    } else {
      ret = tsk_edge_table_update_row(
          &r->tables.edges, edge, r->tables.edges.left[edge], right + 1,
          r->tables.edges.parent[edge], r->tables.edges.child[edge], NULL, 0);
      check_tsk_error(ret);
    }
  }
  for (int i = t->sample_size; i < 2 * t->sample_size - 1; i++) {
    tree_parent = r->tree_to_table_map[i];
    tree_child = r->tree_to_table_map[t->right_child[i]];
    edge = r->right_edge_below[i];
    edge_child = r->tables.edges.child[edge];
    edge_parent = r->tables.edges.parent[edge];
    if (tree_child != edge_child || tree_parent != edge_parent) {
      ret = tsk_edge_table_add_row(&r->tables.edges, right, right + 1,
                                   tree_parent, tree_child, NULL, 0);
      check_tsk_error(ret);
      r->right_edge_below[i] = r->tables.edges.num_rows - 1;
    } else {
      ret = tsk_edge_table_update_row(
          &r->tables.edges, edge, r->tables.edges.left[edge], right + 1,
          r->tables.edges.parent[edge], r->tables.edges.child[edge], NULL, 0);
      check_tsk_error(ret);
    }
  }
  return;
}

void tree_sequence(recorder *r) {
  int ret;
  ret = tsk_table_collection_sort(&r->tables, 0, 0);
  check_tsk_error(ret);
  ret = tsk_table_collection_build_index(&r->tables, 0);
  check_tsk_error(ret);
  ret = tsk_treeseq_init(&r->tree_sequence, &r->tables, 0);
  check_tsk_error(ret);
  return;
}

void free_recorder(recorder *r) {
  tsk_table_collection_free(&r->tables);
  tsk_treeseq_free(&r->tree_sequence);
  free(r->left_edge_below);
  free(r->right_edge_below);
  free(r->tree_to_table_map);
  free(r);
  return;
}
