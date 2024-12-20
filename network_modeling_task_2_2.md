## Task 2



(1) Using the RSiena manual, we write down the formulas for each effect as follows:


* Out-degree $s_{1i}(x) = \sum_j x_{ij}$

* Reciprocity $s_{2i}(x) = \sum_j x_{ij} x_{ji}$

* Transitive reciprocated triplets $s_{3i}(x) = \sum_{j, h} x_{ij} x_{ji} x_{ih} x_{hj}$

* Indegree popularity $s_{4i}(x) = \sum_{j, k, j \neq k} x_{ji} x_{ki}$

* Same covariate $s_{5i}(x, v) = \sum_j x_{ij} \mathbb{I}[v_i = v_j]$

Where $\mathbb{I}$ denotes the indicator function.


(2) We considered it the easiest to write a quick python script that implements each of the effects defined above. We can then run it on the provided graph to obtain complete results (both $f_j$ and $p_j$ values for every possible edge addition), which then makes it easy to answer the provided questions. We share the python script below:

```py
import scipy


# All graph vertices
task_V = {
    'a': 'white',
    'b': 'white',
    'c': 'gray',
    'd': 'gray', 
}

# All directed graph edges
task_E  = [
    ('a', 'b'),
    ('a', 'd'),
    
    ('b', 'c'),
    
    ('c', 'a'),
    
    ('d', 'a'),
    ('d', 'b')
]


# equal to x_{ij} for the given edges
def x(E, i, j):
    return 1 if (i, j) in E else 0

# equal to v_i for the given vertex
def v(V, i):
    return V[i]

# indicator variable
def I(a):
    return 1 if a else 0


# out-degree
def s1i(V, E, i):
    return sum([x(E, i, j) for j in V.keys()])

# reciprocity
def s2i(V, E, i):
    return sum([x(E, i, j) * x(E, j, i) for j in V.keys()])

# transitive reciprocated triplets
def s3i(V, E, i):
    return sum([x(E, i, j) * x(E, j, i) * x(E, i, h) * x(E, h, j) for j in V.keys() for h in V.keys()])

# indegree popularity
def s4i(V, E, i):
    return sum([x(E, j, i) * x(E, k, i) for j in V.keys() for k in V.keys() if j != k])

# same covariate
def s5i(V, E, i):
    return sum([x(E, i, j) * I(V[i] == V[j]) for j in V.keys()])

# the total score of the given vertices and edges
def score(V, E):
    s1 = sum([s1i(V, E, i) for i in V.keys()])
    s2 = sum([s2i(V, E, i) for i in V.keys()])
    s3 = sum([s3i(V, E, i) for i in V.keys()])
    s4 = sum([s4i(V, E, i) for i in V.keys()])
    s5 = sum([s5i(V, E, i) for i in V.keys()])
    
    # compute total score using the provided beta values
    return -1.2 * s1 + 1.5 * s2 + 1 * s3 + 0.5 * s4 + 1.3 * s5





# print the probabilities for all operations for a selected actor i
def step(V, E, i):
    results = []
    # base score (no changes)
    base = score(V, E)
    # go over all possible edges to add / remove
    for j in V.keys():
        if i == j: # source = target -> do nothing action
            results.append(("Do nothing", 0))
        elif not x(E, i, j): # Edge not present in graph - we are adding it
            results.append((f"Adding {i} -> {j}", score(V, (E + [(i, j)])) - base))
        else: # Edge present in graph - we are removing it
            results.append((f"Removing {i} -> {j}", score(V, [e for e in E if e != (i, j)]) - base))
    ops, scores = zip(*results)

    # convert computed scores to probabilities for each action
    next_probs = scipy.special.softmax(scores)

    # print results
    print ("Next step probabilities:" + "".join([f"\n\t* {op}: f={score:.1f}, p={v:.3f}" for op, score, v in zip(ops, scores, next_probs)]))
    
for i in task_V:
    step(task_V, task_E, i)
```





The output is:
```
Next step probabilities:
        * Do nothing: f=0.0, p=0.022
        * Removing a -> b: f=-1.1, p=0.007
        * Adding a -> c: f=3.8, p=0.968
        * Removing a -> d: f=-1.8, p=0.004
Next step probabilities:
        * Adding b -> a: f=8.1, p=0.963
        * Do nothing: f=0.0, p=0.000
        * Removing b -> c: f=1.2, p=0.001
        * Adding b -> d: f=4.8, p=0.036
Next step probabilities:
        * Removing c -> a: f=0.2, p=0.010
        * Adding c -> b: f=4.8, p=0.959
        * Do nothing: f=0.0, p=0.008
        * Adding c -> d: f=1.1, p=0.024
Next step probabilities:
        * Removing d -> a: f=-2.8, p=0.006
        * Removing d -> b: f=0.2, p=0.117
        * Adding d -> c: f=2.1, p=0.782
        * Do nothing: f=0.0, p=0.096
```

We can now easily answer all the subtasks (assuming that each actor was already selected for that particular ministep. To obtain the global value - both sampling the actor and selecting this option - divide percentage values by 4 to account for actor sampling.):
 * Probability that actor c adds a tie to b: $p=95.9\%$
 * Probability that actor b adds a tie to a: $p=96.3\%$
 * Probability that actor a deletes the tie to b: $p=0.7\%$
 * Probability that actor d does not change anything: $p=9.6\%$


