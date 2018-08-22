__precompile__()

"""
CircStats.Datasets includes a number of example circular statistical datasets.
"""
module Datasets

"""Azimuths at which 714 mallards disappeared from sight.

#### Reference
- Matthews, G.V.T. (1961).  "Nonsense" orientation in mallard <em>anas platyrhynchos</em>
  and its relation to experiments on bird navigation.  Ibis, 103a, 211-230.
"""
const Matthews_1961_mallards = (0:20:340, [40, 22, 20, 9, 6, 3, 3, 1, 6, 3, 11, 22, 24, 58, 136, 138, 143, 69])

end # module