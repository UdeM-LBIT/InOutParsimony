# InOutParsimony
Implementation of the InOutParsimony algorithm for solving the Small Gain-Loss Phylogeny problem

The jupyter notebook InOutParsimony.ipynb allows to use the InOutParsimony algorithm to solve the Small Gain-Loss Phylogeny problem. The program takes as input a synteny tree (in newick format), a dictionnary mapping each synteny to a synteny content and costs for gain and loss events. The program outputs a most parsimonious history as a fully labeled tree in Newick format in which gain and loss events are represented as unary node. Each node is labeled by its Name (either the name of the node given in input or Gain/Loss for unary nodes) and its Content (the synteny content at the node for binary nodes and leaves and the content that is gain/loss for gain/loss event).
