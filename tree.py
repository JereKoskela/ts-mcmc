import math
import numpy
import random
import scipy.special


class Tree:
    def __init__(self, sample_size):
        self.left_child = numpy.full(2 * sample_size - 1, -1)
        self.right_child = numpy.full(2 * sample_size - 1, -1)
        self.parent = numpy.full(2 * sample_size - 1, -1)
        self.time = numpy.zeros(2 * sample_size - 1)
        self.sample_size = sample_size

        t = 0
        active_lineages = list(range(sample_size))
        n = len(active_lineages)
        next_parent = sample_size
        while n > 1:
            t += random.expovariate(lambd=scipy.special.binom(n, 2))
            [l_child, r_child] = numpy.random.choice(
                active_lineages, size=2, replace=False
            )
            self.parent[l_child] = next_parent
            self.parent[r_child] = next_parent
            self.left_child[next_parent] = l_child
            self.right_child[next_parent] = r_child
            self.time[next_parent] = t
            active_lineages.remove(l_child)
            active_lineages.remove(r_child)
            active_lineages.append(next_parent)
            next_parent += 1
            n -= 1

    def sample_reattach_time(self, child, new_sib):
        sib = self.sibling(child)
        new_sib_parent = self.parent[new_sib]
        new_sib_grandparent = self.grandparent(new_sib)
        if sib == new_sib:
            if new_sib_grandparent == -1:
                ret = self.time[new_sib] + random.expovariate(lambd=1)
            else:
                ret = numpy.random.uniform(
                    self.time[new_sib], self.time[new_sib_grandparent]
                )
        else:
            if new_sib_parent == -1:
                ret = self.time[new_sib] + random.expovariate(lambd=1)
            else:
                ret = numpy.random.uniform(
                    self.time[new_sib], self.time[new_sib_parent]
                )
        return ret

    def log_reattach_density(self, child, new_sib, new_time):
        sib = self.sibling(child)
        new_sib_parent = self.parent[new_sib]
        new_sib_grandparent = self.grandparent(new_sib)
        if sib == new_sib:
            if new_sib_grandparent == -1:
                ret = self.time[new_sib] - new_time
            else:
                ret = -math.log(self.time[new_sib_grandparent] - self.time[new_sib])
        else:
            if new_sib_parent == -1:
                ret = self.time[new_sib] - new_time
            else:
                ret = -math.log(self.time[new_sib_parent] - self.time[new_sib])
        return ret

    def resample_times(self):
        ret = 0
        lb = 0
        [sorted_times, ind] = numpy.unique(self.time, return_index=True)
        for i in range(1, len(sorted_times)):
            rate = scipy.special.binom(self.sample_size - i + 1, 2)
            index = ind[i]
            self.time[index] = lb + random.expovariate(lambd=rate)
            ret -= rate * (self.time[index] - lb)
            lb = self.time[index]
        return ret

    def log_resample_times_density(self, t):
        ret = 0
        lb = 0
        [sorted_times, ind] = numpy.unique(t, return_index=True)
        for i in range(1, len(sorted_times)):
            rate = scipy.special.binom(self.sample_size - i + 1, 2)
            index = ind[i]
            ret -= rate * (t[index] - lb)
            lb = t[index]
        return ret

    def sample_leaf(self):
        node = numpy.random.choice(numpy.arange(self.sample_size))
        return node

    def sample_node(self):
        node = numpy.random.choice(len(self.parent))
        return node

    def replace_child(self, node, child, new_child):
        if self.left_child[node] == child:
            self.left_child[node] = new_child
        else:
            self.right_child[node] = new_child

    def detach_reattach(self, child, new_sib, new_time):
        sib = self.sibling(child)
        child_parent = self.parent[child]
        sib_grandparent = self.grandparent(child)
        new_sib_parent = self.parent[new_sib]
        self.time[child_parent] = new_time
        if sib != new_sib:
            self.parent[new_sib] = child_parent
            self.parent[sib] = sib_grandparent
            self.parent[child_parent] = new_sib_parent
            self.replace_child(child_parent, sib, new_sib)
            if sib_grandparent != -1:
                self.replace_child(sib_grandparent, child_parent, sib)
            if new_sib_parent != -1:
                self.replace_child(new_sib_parent, new_sib, child_parent)

    def sibling(self, node):
        ret = -1
        p = self.parent[node]
        if p != -1:
            ret = self.left_child[p]
            if ret == node:
                ret = self.right_child[p]
        return ret

    def grandparent(self, node):
        ret = -1
        if self.parent[node] != -1:
            ret = self.parent[self.parent[node]]
        return ret

    def log_likelihood(self):
        sorted_times = numpy.unique(self.time)
        ret = 0
        for i in range(self.sample_size - 1):
            ret -= scipy.special.binom(self.sample_size - i, 2) * (
                sorted_times[i + 1] - sorted_times[i]
            )
        return ret
