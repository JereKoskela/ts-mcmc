import math
import numpy
import random

from recorder import Recorder
from tree import Tree


def kingman_mcmc(tree, recorder):
    acceptance_prob_spr = 0
    acceptance_prob_times = 0
    log_likelihood = tree.log_likelihood()
    steps = int(recorder.tables.sequence_length) - 1

    for i in range(steps):
        child = tree.sample_leaf()
        parent = tree.parent[child]
        sib = tree.sibling(child)
        new_sib = sib
        if tree.sample_size > 2:
            while new_sib in [child, parent, sib]:
                new_sib = tree.sample_node()
        new_time = tree.sample_reattach_time(child, new_sib)
        alpha = -tree.log_reattach_density(child, new_sib, new_time)
        old_time = tree.time[parent]
        tree.detach_reattach(child, new_sib, new_time)
        alpha += tree.log_reattach_density(child, sib, old_time)
        proposal_log_likelihood = tree.log_likelihood()
        alpha += proposal_log_likelihood - log_likelihood
        if math.log(random.random()) < alpha:
            log_likelihood = proposal_log_likelihood
            acceptance_prob_spr += 1 / steps
        else:
            tree.detach_reattach(child, sib, old_time)

        old_times = tree.time.copy()
        alpha = -tree.resample_times()
        proposal_log_likelihood = tree.log_likelihood()
        alpha += tree.log_resample_times_density(old_times)
        alpha += proposal_log_likelihood - log_likelihood
        if math.log(random.random()) < alpha:
            log_likelihood = proposal_log_likelihood
            acceptance_prob_times += 1 / steps
        else:
            tree.time = old_times.copy()

        recorder.append_tree(tree)

    return [acceptance_prob_spr, acceptance_prob_times]
