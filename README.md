# Simplicial Complex in Communication Networks
---
- The communication networks used in the paper are located in the data folder. However, due to their large size, the phone calls communication network (https://www.science.org/doi/10.1126/science.1177170) and the University of Kiel email network (http://www.theo-physik.uni-kiel.de/%7Eebel/email-net/email_net.html) are not included in the GitHub repository.
- The random networks that have been generated are located in the Output folder.
- `enronSimplices.py` is used to load the networks and extract S*- and T*-simplices from them.
- `para2para.py` performs symbolic regression on the distribution of simplicial complexes and generates visualizations of the estimated parameters.
- `randomNetwork.py` creates a random power-law network based on the specified n, m, p settings.
- `shuffledEnron.py` analyzes the distribution of simplices after randomizing the original Enron network.
- `SimplicesExtract.py` extracts qualitative simplices from the Enron network as discussed in the paper.