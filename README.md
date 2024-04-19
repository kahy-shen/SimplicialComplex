# Simplicial Complex in Communication Networks
---
- The communication networks employed in the paper are stored in the data folder.
- The random networks that have been generated are located in the Output folder.
- `enronSimplices.py` is used to load the networks and extract S*- and T*-simplices from them.
- `para2para.py` performs symbolic regression on the distribution of simplicial complexes and generates visualizations of the estimated parameters.
- `randomNetwork.py` creates a random power-law network based on the specified n, m, p settings.
- `shuffledEnron.py` analyzes the distribution of simplices after randomizing the original Enron network.
- `SimplicesExtract.py` extracts qualitative simplices from the Enron network as discussed in the paper.