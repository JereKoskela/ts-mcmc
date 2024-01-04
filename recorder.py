import math
import numpy
import random
import tskit

from tree import Tree


class Recorder:
    def __init__(self, tree, sample_size, seq_length):
        self.tables = tskit.TableCollection(sequence_length=seq_length)
        self.node_table = self.tables.nodes
        self.edge_table = self.tables.edges
        num_nodes = len(tree.parent)
        self.active_edges = list(range(num_nodes - 1))
        self.tree_to_table_map = list(range(num_nodes))
        for i in self.tree_to_table_map:
            if i < sample_size:
                node = self.node_table.add_row(time=tree.time[i], flags=1)
            else:
                node = self.node_table.add_row(time=tree.time[i])
            if tree.left_child[i] != -1:
                self.edge_table.add_row(
                    left=0, right=1, child=tree.left_child[i], parent=i
                )
                self.edge_table.add_row(
                    left=0, right=1, child=tree.right_child[i], parent=i
                )
        self.tables.sort()
        ts = self.tables.tree_sequence()

    def store_spr(self, tree, child, sib):
        r = self.edge_table[self.active_edges[0]].right
        # tree has already undergone the MCMC update so parent etc. indices
        # need to be retrofitted.
        parent = tree.parent[child]
        new_sib = tree.sibling(child)
        if new_sib == sib:
            new_sib_parent = parent
            sib_grandparent = tree.grandparent(sib)
        else:
            new_sib_parent = tree.grandparent(new_sib)
            sib_grandparent = tree.parent[sib]
        ts_parent = self.tree_to_table_map[parent]

        ts_new_parent = self.node_table.add_row(time=tree.time[parent])
        self.tree_to_table_map[parent] = ts_new_parent
        ts_child = self.tree_to_table_map[child]
        next_edge = self.edge_table.add_row(
            left=r, right=r + 1, child=ts_child, parent=ts_new_parent
        )
        ts_new_sib = self.tree_to_table_map[new_sib]
        self.edge_table.add_row(
            left=r, right=r + 1, child=ts_new_sib, parent=ts_new_parent
        )
        new_edges = 2
        ts_sib_grandparent = self.tree_to_table_map[sib_grandparent]
        ts_sib = self.tree_to_table_map[sib]
        ts_new_sib_parent = self.tree_to_table_map[new_sib_parent]
        if new_sib == sib:
            if sib_grandparent != -1:
                self.edge_table.add_row(
                    left=r, right=r + 1, child=ts_new_parent, parent=ts_sib_grandparent
                )
                new_edges += 1
        else:
            new_edges += 1
            if sib_grandparent == -1:
                self.edge_table.add_row(
                    left=r, right=r + 1, child=ts_new_parent, parent=ts_new_sib_parent
                )
            elif new_sib_parent == -1:
                self.edge_table.add_row(
                    left=r, right=r + 1, child=ts_sib, parent=ts_sib_grandparent
                )
            else:
                self.edge_table.add_row(
                    left=r, right=r + 1, child=ts_sib, parent=ts_sib_grandparent
                )
                self.edge_table.add_row(
                    left=r, right=r + 1, child=ts_new_parent, parent=ts_new_sib_parent
                )
                new_edges += 1

        for edge in reversed(self.active_edges):
            if self.edge_table[edge].child in [
                ts_child,
                ts_parent,
                ts_sib,
                ts_new_sib,
            ]:
                self.active_edges.remove(edge)
            else:
                self.edge_table[edge] = self.edge_table[edge].replace(right=r + 1)
        for i in range(new_edges):
            self.active_edges.append(next_edge)
            next_edge += 1

    def increment_site(self):
        r = self.edge_table[self.active_edges[0]].right
        for edge in self.active_edges:
            self.edge_table[edge] = self.edge_table[edge].replace(right=r + 1)

    def tree_sequence(self):
        self.tables.sort()
        ts = self.tables.tree_sequence()
        return ts
