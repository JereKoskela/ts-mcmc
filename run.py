from mcmc import kingman_mcmc
from recorder import Recorder
from tree import Tree

sample_size = 4
steps = 10
tree = Tree(sample_size)
recorder = Recorder(tree, sample_size, steps)
acceptance_prob = kingman_mcmc(tree, recorder)
ts = recorder.tree_sequence()
print(ts.draw_text())
print("acceptance probability = ", acceptance_prob)
