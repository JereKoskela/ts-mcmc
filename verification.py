import numpy
import matplotlib.pyplot as plt
import tskit

from mcmc import kingman_mcmc
from recorder import Recorder
from tree import Tree


def verify_kingman(samples, steps):
    len_samples = len(samples)
    mean_tmrca = numpy.zeros(len_samples)
    var_tmrca = numpy.zeros(len_samples)
    mean_branch_length = numpy.zeros(len_samples)
    var_branch_length = numpy.zeros(len_samples)
    acceptance_probs = numpy.zeros(len_samples)
    ind = 0
    for i in samples:
        tree = Tree(i)
        recorder = Recorder(tree, i, steps)
        acceptance_prob = kingman_mcmc(tree, recorder)
        ts = recorder.tree_sequence()
        node_table = recorder.node_table
        for tree in ts.trees():
            mean_tmrca[ind] += (
                node_table[tree.root].time
                * (tree.interval.right - tree.interval.left)
                / steps
            )
            var_tmrca[ind] += (
                node_table[tree.root].time ** 2
                * (tree.interval.right - tree.interval.left)
                / steps
            )
            mean_branch_length[ind] += (
                tree.total_branch_length
                * (tree.interval.right - tree.interval.left)
                / steps
            )
            var_branch_length[ind] += (
                tree.total_branch_length**2
                * (tree.interval.right - tree.interval.left)
                / steps
            )
        var_tmrca[ind] -= mean_tmrca[ind] ** 2
        var_branch_length[ind] -= mean_branch_length[ind] ** 2
        acceptance_probs[ind] = acceptance_prob
        ind += 1
    exact_mean_tmrca = numpy.array([2 * (1 - 1 / n) for n in samples])
    exact_var_tmrca = numpy.cumsum([(2 / (n * (n - 1))) ** 2 for n in samples])
    exact_mean_branch_length = numpy.cumsum([2 / (n - 1) for n in samples])
    exact_var_branch_length = numpy.cumsum([4 / (n - 1) ** 2 for n in samples])

    plt.plot(samples, mean_tmrca, label="Observed E[TMRCA]")
    plt.plot(
        samples,
        exact_mean_tmrca,
        label="Expected E[TMRCA]",
    )
    plt.plot(samples, var_tmrca, label="Observed Var(TMRCA)")
    plt.plot(
        samples,
        exact_var_tmrca,
        label="Expected Var(TMRCA)",
    )
    plt.legend()
    plt.savefig("tmrca.png")
    plt.close("all")

    plt.plot(
        samples,
        mean_branch_length,
        label="Observed E[Branch length]",
    )
    plt.plot(
        samples,
        exact_mean_branch_length,
        label="Expected E[Branch length]",
    )
    plt.plot(
        samples,
        var_branch_length,
        label="Observed Var(Branch length)",
    )
    plt.plot(
        samples,
        exact_var_branch_length,
        label="Expected Var(Branch length)",
    )
    plt.legend()
    plt.savefig("branch_length.png")
    plt.close("all")

    print(acceptance_probs)


min_sample = 2
max_sample = 10
samples = range(min_sample, max_sample + 1)
steps = int(1e6)

verify_kingman(samples, steps)
