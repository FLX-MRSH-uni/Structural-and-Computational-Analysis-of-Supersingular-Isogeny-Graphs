GL_graph_analysis.sage

This document is the accompanying code for the year 4 MSCi Mathematics dissertation "Structural and Computational Analysis of
Supersingular Isogeny Graphs". Analysis of some experimental results using this code document can be found in said dissertation.
The reader is encouraged to read the code before using it in experiments. Please note this is a Sage file and not a Python file. 
Here is a guide to the three important functions:


.range_pL_eigenvalue_data(n, m, L, highlight_p = None, plot_lambda = True, plot_quantity = True) 

For all primes p ≡ 1 (mod 12) between n and m, computes the eigenvalues of the adjacency matrix of each GL graph, and catagorises  
them into (all) non-trivial eigenvalues, and eigenvalues with absolute value greater than the un-Ramanujan bound  2√(d − 1).
plot_lambda plots the largest such absolute value and plot_quantity plots the quantity of such eigenvalues.
Returns eigenvalue data.


.compare_L_dist(p, L1, L2, bin_width=0.0001, plot=True)

Computes the shortest path distributions, by means of time complexity (estimate of the time taken to compute the isogenies in the path),
with specific p for multiple GL with L1 ≤ L ≤ L2 in lexicographical order. bin_width is the range between time complexities 
represented by each data point. The distributions are plotted on one another for comparison.
Returns largest recorded difference in diamter. 


.compute_and_plot_MitM_data(p, L1, L2, fixed_sens=None, fix_pairs=False, fixed_no_pairs=8, plot=True)

Computes a mean a standard deviation for estimates of the time complexity of a Meet in the Middle attack on GL graphs for specific p 
with L1 ≤ L ≤ L2 in lexicographical order and plots the data against un-Ramanujanness. Pairs of vertices used are either fixed for each
L (fix_pairs), with a fixed number used, or are "average pairs" for each GL graph by means of proximity to the mean shortest path of 
each graph. The proximity to the mean (sensitivity) can be fixed if desired. 
Returns un-Ramanujanness and statistical MitM data for each graph, as well as the dominant value of L and the sensitivity used. 
