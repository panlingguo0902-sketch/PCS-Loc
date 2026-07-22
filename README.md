# PCS-Loc
rumor source
propagation consistency-based
rumor source localization approach

PCS-Loc algorithm is designed for estimating the rumor source node in social networks. The PCS-Loc algorithm requires the following installation setup.

Python 3.10 are above
NDlib - Network Diffusion Library
The Python code (PCS-Loc) requires the following packages to be imported

import networkx as nx
import csv
import random
import ndlib.models.ModelConfig as mc
import ndlib.models.epidemics as ep
import numpy as np
import timeit
from collections import deque
This code requires an Input file (.csv format) as shown in the following format: (It has two columns such as Source and Target)

Source,Target

1,2

2,1

1,3

3,1

2,3

3,2

The original source can be obtained from the simulation part (IC model).

The output of the ROSE algorithm code is the estimated source, i.e., any node belongs to the graph (input_file.csv).

Once the estimated source and the original source are available, we can generate Distance error, Execution time, Accuracy, Candidate sources, etc.
