# PCS-Loc
rumor source
ATechnique for Identifying Rumor Origins
Within Communication Consistency
Limitations

The PCS-Loc algorithm is proposed to identify the rumor source node in social networks. To run the PCS-Loc algorithm, the corresponding environment configuration and library installation are required as follows.
Prerequisites:
Python 3.10 or higher
NDlib (Network Diffusion Library)
The Python implementation of PCS-Loc requires the following modules to be imported:
import networkx as nx
import csv
import random
import ndlib.models.ModelConfig as mc
import ndlib.models.epidemics as ep
import numpy as np
import timeit
from collections import deque
The program takes a CSV file as input. This file consists of two columns named Source and Target, which describe the edges of the network. An example of the input format is displayed below:
Source,Target
1,2
2,1
1,3
3,1
2,3
3,2
The ground-truth rumor source is generated from simulations based on the IC model.
The PCS-Loc algorithm outputs an estimated rumor source node, which must be a valid node contained in the graph defined by the input CSV file.
With both the estimated source and ground-truth source obtained, multiple evaluation metrics can be computed, including distance error, execution time, detection accuracy, candidate source set and other results.
