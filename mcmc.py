import math
import numpy
import random

from recorder import Recorder
from tree import Tree


def kingman_mcmc(tree, recorder):
    acceptance_prob = 0
    log_likelihood = tree.log_likelihood()
    steps = int(recorder.tables.sequence_length) - 1

    for i in range(steps):
        child = tree.sample_leaf()
        parent = tree.parent[child]
        sib = tree.sibling(child)

        new_sib = tree.sample_node()
        while new_sib in [child, parent]:
            new_sib = tree.sample_node()

        new_time = tree.sample_proposal(child, new_sib)
        alpha = -tree.log_proposal_density(child, new_sib, new_time)
        old_time = tree.time[parent]
        tree.detach_reattach(child, new_sib, new_time)
        alpha += tree.log_proposal_density(child, sib, old_time)
        proposal_log_likelihood = tree.log_likelihood()
        alpha += proposal_log_likelihood - log_likelihood

        if math.log(random.random()) < alpha:
            recorder.store_spr(tree, child, sib)
            log_likelihood = proposal_log_likelihood
            acceptance_prob += 1 / steps
        else:
            tree.detach_reattach(child, sib, old_time)
            recorder.increment_site()
    return acceptance_prob
