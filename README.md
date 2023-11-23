# Fast Multilayer Core Decomposition and Indexing

The multilayer (ML) graph model provides a robust representation of multi-sourced relationships among real-world entities, laying a solid foundation for reliable knowledge discovery. ML core decomposition is a fundamental analytical tool for ML graphs. It offers valuable insights into the dense structures in ML graphs and forms the basis for many complex analysis tasks. However, existing ML core decomposition algorithms face performance issues due to unavoidably unnecessary computations and are inherently serial, unable to fully leverage multi-core processors.
To address this, we reformulate the search space of this problem with a tree-shaped structure called MLC-tree. Based on it, we present an efficient serial ML core decomposition algorithm that improves the time complexity over existing solutions and the first parallel framework for this problem by exploiting the path-decomposition of the MLC-tree. Two practical optimizations are introduced to further boost the parallel efficiency.
To facilitate applications built upon ML cores, we construct a compact storage and index structure for ML cores based on the MLC-tree. The utility of this index is showcased through two applications: ML core search and a novel weighted densest subgraph discovery problem.


<!--Extensive experiments on $9$ real-world ML graphs show that our MLC-tree-based ML core decomposition algorithm achieves a speedup of up to $128\times$ over existing baselines and the parallel approach attains an additional speedup of up to $30.6\times$ using $40$ cores. Moreover, the MLC-tree index proves to efficiently support the studied applications.
## Execution -->

### Usage

```commandline
MlcDec mode [options]
```
#### Mode:

    bt      Build the augmented MLC-tree (include the process of computing the ML core decomposition)
    bht     Build the hashtable holding the ML core decomposition
    bdt     Build the edge-difference-augmented MLC-tree
    
    mlcs    Peeling-based ML core search
    mlcs_t  MLC-tree-based ML core search
    mlcs_h  Hashtable-based ML core search
    
    ds      Weighted densest subgraph discovery
    dts     Weighted densest subgraph discovery based on the edge-difference-augmented MLC-tree

	grk     Generate random \mathbf{k} vectors
	
#### Options:
    -g      Graph name
    -o      Layer Ordering
    -o      Output path
    -ntc	Number of testcases
    -kvf    File of queried coreness vectors
    -mlctf  MLC-tree file
    -mlctf  Hashtable file
    
    -m      Methods of running MLC dec
    -t      Number of threads
    -l      Levels of applying the core-level-parallel startup
    -alpha  Parameter alpha in path merging
    -f      Flush the MLC-tree/hashtable to file
    
    -beta   Tradeoff in the weighted densest subgraph
    -w      Vector of weights

#### Graph ordering options:
    d       Non-desceasing order of degeneracy
    d_r     Non-increasing order of degeneracy
    c       Non-desceasing order of coreness
    c_r     Non-increasing order of degeneracy
    rd      Random order
    


#### Methods of MLC Dec:
    s/serial       Serail
    ip             Basic framework with path-level parallelism
    ips            ip + core-level-parallel startup
    ipm            ip + path merging
    opt            ip + all optimizations

#### Examples: 
    
    'MlcDec bt -g dblp -o d'
    'MlcDec bt -g dblp -o -m ip -t 8'
    'MlcDec bt -g dblp -o -m opt -alpha 0.1 -l 2 -t 8'
    
    'MlcDec mlcs -kvf vector.file'
    'MlcDec mlcs_t -kvf vector.file -mlctf mlct.file'
    
    'MlcDec ds -g dblp -alpha 0.1 -l 2 -beta 2 -w [1,2,3,4,5,6,7,8,9,10] -t 8'
    'MlcDec dts -mlctf ea-mlct.file -beta 2 -w [1,2,3,4,5,6,7,8,9,10]'
    
    