"""
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

"""

# =================================================================================================================================================

import os
import json
import time
import math
import statistics

# =================================================================================================================================================

# This is the code provided in "Cycles and Cuts in Supersingular L-Isogeny Graphs" [1] which we will use for our experiments. 

def get_L_graph(p, L) -> DiGraph:
    """
    Computes the L-isogeny graph over Fp2 for the given list of L's.
    - Each edge is of the form (start_vertex, end_vertex, "{degree},{i}) where i is a counter that is useful to distinguish multi-edges
    """
    Fp2 = GF(p ** 2)
    j = supersingular_j(Fp2)
    E = EllipticCurve_from_j(j)
    edges = []
    for ell in L:
        G_ell = E.isogeny_ell_graph(ell, directed = True, label_by_j = True)      
        seen = []
        for edge in G_ell.edges():
            count  = seen.count(edge)
            edges.append((edge[0], edge[1], f"{ell},{count}"))               
            seen.append(edge)
    return DiGraph(edges, loops = True, multiedges = True)    

# =================================================================================================================================================

# This is the class which contains all the functions used in the experiments presented in the paper,
# along with all the helpers which build up to the functions

class LIsogenyGraphs:

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------------------------------------------------------------------

    # First the initialisation of the data files 

    def __init__(self, data_file, data_file2, data_file3, data_file4, data_file5, data_file6):
        """
        When initialising the class, keep in mind the order of the files:
        - data_file will contain the average time complexities of specific ℓ isogenies from each vertex in the graph
        - data_file2 will contain the distance between each pair of vertices (in terms of time complexity)
        - data_file3 will contain all the non-trivial eigenvalue data for each graph 
        - data_file4 will contain the experimental meet in the middle data 
        - data_file5 will contain the analysis of the raw data from the other files, i.e. the data used for the plots
        - data_file6 will contain the "average pairs" used for MitM 
        """

        self.data_file = data_file
        self.data_file2 = data_file2
        self.data_file3 = data_file3
        self.data_file4 = data_file4
        self.data_file5 = data_file5
        self.data_file6 = data_file6

        self.path_times = self._load_json(self.data_file)
        self.shortest_path_times = self._load_json(self.data_file2)
        self.eigenvalue_data = self._load_json(self.data_file3)
        self.MitM_data = self._load_json(self.data_file4)
        self.analysis = self._load_json(self.data_file5)
        self.average_pairs = self._load_json(self.data_file6)

    # This is a helper for opening the data files 

    def _load_json(self, path):
        if not os.path.exists(path):
            return {}
        try:
            with open(path, "r") as f:
                return json.load(f)

        # be careful when cancelling any computation as you might truncate the JSON file  
        # and make it unreadable to the JSON module, 
        # if you do cancel, look at the files and edit accordingly as to not loose your data

        except json.JSONDecodeError:
            print(f"WARNING: {path} is not valid JSON (likely truncated). Reinitialising to empty dict.")
            return {}

    # This is a helper for saving the data files

    def _save_json(self, path, obj):

        obj = self.to_jsonable(obj)

        if isinstance(obj, dict):
            try:
                obj = dict(sorted(obj.items(), key=lambda kv: int(kv[0])))
            except Exception:
                pass

        with open(path, "w") as f:
            json.dump(obj, f, indent=2)

    # This is a helper which sorts data by the choice of L in lexicographical order

    def _sort_L_keys(self, d):
        return dict(sorted(d.items(), key=lambda kv: json.loads(kv[0])))

    # This is a helper hich converts the choice of L into a string readable in a JSON file

    def _L_key(self, L):
        L_int = sorted(int(x) for x in L)
        return json.dumps(L_int, separators=(",", ":"))

    # This converts various data types into a string format for a JSON file

    def to_jsonable(self, x):
        if isinstance(x, dict):
            return {str(k): self.to_jsonable(v) for k, v in x.items()}

        if isinstance(x, (list, tuple)):
            return [self.to_jsonable(v) for v in x]

        if isinstance(x, float):
            return float(x)

        if x.__class__.__name__ == "Integer":
            return int(x)

        if x.__class__.__name__ in ("RealNumber", "RealLiteral"):
            return float(x)

        try:
            if hasattr(x, "n"):
                return float(x)
        except Exception:
            pass

        try:
            if hasattr(x, "is_integer") and x.is_integer():
                return int(x)
        except Exception:
            pass

        return x

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Experiment 1: Eigenvalues
# -------------------------------------------------------------------------------------------------------------------------------------------------

    # First, eigenvalue data for a single graph

    def pL_eigenvalue_data(self, p, L):
        """
        For choices of p with p ≡ 1 (mod 12) computes the non-trivial 
        eigenvalues of the adjacency matrix of the GL graph, and formats 
        them into all eigenvalues and those with abs above the threshold.
        Returns:
        - Non-trivial eigenvalues
        - Above threshold eigenvalues
        """

        # Ensure p ≡ 1 (mod 12)

        if not p % 12 == 1:
            print(f"WARNING: {p} is not 1(mod 12).")
            return None

        # form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)

        # load or create the specific area in the eigenvalue data file

        self.eigenvalue_data.setdefault(p_key, {})
        self.eigenvalue_data[p_key].setdefault("L_specific", {})
        self.eigenvalue_data[p_key]["L_specific"].setdefault("eigenvalue_data", {})

        eig_data = self.eigenvalue_data[p_key]["L_specific"]["eigenvalue_data"]

        # check if the data already exists in our data file and if it does we return it

        if L_key in eig_data:
            return eig_data[L_key]

        # otherwise build the graph, compute the adjacency matrix, and compute the eigenvalues

        G = get_L_graph(p, L)
        A = G.adjacency_matrix()
        eigenvals = A.eigenvalues()
        eigenvals_RR = [RR(ev) for ev in eigenvals]

        # find which eigenvalue has the largest absolute value (the trivial eigenvalue)

        max_abs = 0
        max_ind = 0
        for i in range(len(eigenvals_RR)):
            val = eigenvals_RR[i]
            if abs(val) > max_abs:
                max_abs = abs(val)
                max_ind = i

        # remove the trivial eigenvalue and sort the rest of them by absolute value

        non_trivial_eigs = [eigenvals_RR[i] for i in range(len(eigenvals_RR)) if i != max_ind]
        non_trivial_eigs_sorted = sorted(non_trivial_eigs, key=lambda x: abs(x), reverse=True)
        non_trivial_eigs_sorted = [float(x) for x in non_trivial_eigs_sorted]

        # compute the un-Ramanujan threshold and determine which eigenvalues are greater than this threshold

        d = sum(L)+len(L)
        threshold = 2 * sqrt(d - 1)
        violating_eigs_sorted = [lam for lam in non_trivial_eigs_sorted if abs(lam) > threshold]
        violating_eigs_sorted = [float(x) for x in violating_eigs_sorted]

        # format the data

        entry = {
            "nontrivial": non_trivial_eigs_sorted,
            "above_threshold": violating_eigs_sorted
        }

        # save the data to the data file and return it

        eig_data[L_key] = entry
        self.eigenvalue_data[p_key]["L_specific"]["eigenvalue_data"] = dict(
            sorted(
                self.eigenvalue_data[p_key]["L_specific"]["eigenvalue_data"].items(),
                key=lambda kv: json.loads(kv[0])
            )
        )
        self._save_json(self.data_file3, self.eigenvalue_data)
        return entry

    # Now, eigenvalue data for a range of primes so we can visualise any potential trend

    def range_pL_eigenvalue_data(self, n, m, L, highlight_p = None, plot_lambda = True, plot_quantity = True):
        """
        For primes p ≡ 1 (mod 12) between n and m and a chosen L, compute the non-trivial eignvalue data for each GL graph and 
        - if plot_lambda = True: plot the largest such absolute value for each graph with the un-Ramanujan threshold visible,
        - if plot_quantity = True: plot the number of eigenvalues with abs greater than the threshold,
        - if a list of primes highlight_p is given: will highlight those specific primes on each plot.
        Returns:
        - Dictionary containing non-trivial eigenvalues of each p,
        - Dictionary containing every eigenvalue with abs above threshold for each p.
        """
        
        # Ensure p ≡ 1 (mod 12)

        primes = [p for p in prime_range(n+1, m) if p % 12 == 1]

        # Compute the threshold (this is only dependent on L)

        d = sum(L)+len(L)
        threshold = 2 * sqrt(d - 1)

        # initialise the relevant dictionaries and lists

        non_trivial_res = {}
        violating_res = {}
        counts = []
        largest = []

        # this for loop fills up the above with the relevant data using the previous function,
        # if the data already exists the previous function will find and return it, 
        # otherise it will calculate it

        for p in primes:
            p_data = self.pL_eigenvalue_data(p, L)
            non_trivial_res[p] = p_data["nontrivial"]
            violating_res[p] = p_data["above_threshold"]
            counts.append((p, len(p_data["above_threshold"])))
            if len(p_data["nontrivial"]) > 0:
                largest.append((p, abs(p_data["nontrivial"][0])))

        # if there is a list of primes to highlight, seperate them from the others

        highlight_set = set(int(q) for q in (highlight_p or []))
        largest_hi = [pt for pt in largest if int(pt[0]) in highlight_set]
        largest_lo = [pt for pt in largest if int(pt[0]) not in highlight_set]
        counts_hi  = [pt for pt in counts  if int(pt[0]) in highlight_set]
        counts_lo  = [pt for pt in counts  if int(pt[0]) not in highlight_set]

        # convert L into a string readable for plotting purposes

        L_id = "".join(str(x) for x in L)

        # plot the data 

        if plot_lambda and primes:

            # first plot the non-highlighted

            K = list_plot(
                largest_lo,
                plotjoined=False,
                marker='o'
            )
            
            # then the highlighted

            if largest_hi:
                K += list_plot(largest_hi, plotjoined=False, marker='o', color='brown')

            xmin = min(primes)
            xmax = max(primes)
            line_plot = line([(xmin, threshold), (xmax, threshold)], color='red', thickness=1.5, zorder=0)

            # then add the threshold

            K_with_line = line_plot + K

            # format the plot and add the labels

            figsize = (6.3, 5.2)
            fig = K_with_line.matplotlib(figsize=figsize)
            ax = fig.axes[0]
            ax.set_xlabel(r'$p$', fontsize=12)
            ax.set_ylabel(r'$max(|\lambda|)$', fontsize=12)

            # save the figure

            filename = f"largest_lambda_L_isogeny_plot_{min(primes)}_{max(primes)}_{L_id}.png"
            fig.savefig(filename)
            print(f"Saved plot to {filename}")

        # plot as above

        if plot_quantity and primes:
            P = list_plot(
                counts_lo,
                plotjoined=False,
                marker='o'
            )

            if counts_hi:
                P += list_plot(counts_hi, plotjoined=False, marker='o', color='brown')

            figsize = (6.3, 5.2)
            fig = P.matplotlib(figsize=figsize)
            ax = fig.axes[0]
            ax.set_xlabel(r'$p$', fontsize=12)
            ax.set_ylabel(r'$\#\{\lambda: |\lambda| > 2\sqrt{d-1}\}$')
            filename = f"no_large_lambda_L_isogeny_plot_{min(primes)}_{max(primes)}_{L_id}.png"
            fig.savefig(filename)
            print(f"Saved plot to {filename}")

        # return all eigenvalues and all eigenvalues with abs above the threshold

        return non_trivial_res, violating_res

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Experiment 2: Diameter and Shortest Paths
# -------------------------------------------------------------------------------------------------------------------------------------------------

    # Determine the edge weights in each GL graph by timing computation

    def isogeny_costs(self, p, L):
        """
        Determine the average time complexity of each isogeny degree ℓ for each vertex in GL 
        by timing total computation for each ℓ and dividing by the total number of outgoing ℓ-isogenies.
        Returns:
        - Dictionary containing the average time complexity of each isogeny degree ℓ for each vertex in GL.
        """

        # Form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)

        # load or create the specific area in the isogeny costs data file

        self.path_times.setdefault(p_key, {})
        pL_data = self.path_times[p_key]

        # if the data already exists, return it 

        if L_key in pL_data:
            return pL_data[L_key]

        # this time initialise the field as we need it for the isogeny computations

        Fp2 = GF(p^2)

        # the ℓ-isogenies in a GL graph are the same as the ones in a Gℓ graph,
        # to make the function less convoluted, we can simply treat each ℓ as seperate cases,
        # thus timing the isogenies in each Gℓ graph, and compiling the data for the GL case

        for ell in L:

            # if the data already exists, don't recompute

            ell_key = "["+str(int(ell))+"]"
            if ell_key in pL_data:
                continue

            # create the area in the data file

            self.path_times[p_key].setdefault(ell_key, {})

            # make the Gℓ graph

            G_ell = get_L_graph(p, [ell])

            # now for each vertex in Gℓ,

            for v in G_ell.vertices():

                # make sure it has neighbours (which it should),

                neighbours = [w for w in G_ell.neighbor_out_iterator(v)]
                if not neighbours:
                    continue

                # convert the vertex into an elliptic curve in the isomorphism class of the j-invariant,

                j_val = Fp2(v)
                E_v = EllipticCurve_from_j(j_val)

                # use the time module to time how long it takes to compute every ℓ-isogeny from the elliptic curve,
                # here we use the sage function .isogenies_prime_degree() to find the list of such isogenies,
                # see [2]

                t0 = time.time()
                isos = E_v.isogenies_prime_degree(ell)
                t1 = time.time()

                dt = t1 - t0
                if not isos:
                    continue

                # find the average for the vertex

                cost_per_isogeny = dt / len(isos)

                # save the data for the Gℓ graph

                pL_data[ell_key][str(v)] = [(ell, cost_per_isogeny)]
                self._save_json(self.data_file, self.path_times)
                
        # now compile the data from each Gℓ graph into data for the GL graph 

        avg_times = {}
        for ell in L:
            ell_key = "["+str(int(ell))+"]"

            for v in pL_data[ell_key]:
                if v in avg_times.keys():
                    avg_times[v] += [(ell, pL_data[ell_key][v][0][1])]
                else:
                    avg_times[v] = [(ell, pL_data[ell_key][v][0][1])]                  

        # make sure the area in the data file exists

        if not L_key in pL_data:
            pL_data[L_key] = avg_times

        # save and return the data

        self.path_times[p_key] = self._sort_L_keys(self.path_times[p_key])
        self._save_json(self.data_file, self.path_times)
        return avg_times                               

    # This helper determines the degree of each edge (value of ℓ)

    def edge_degree(self, edge):
        ell, _ = edge[2].split(",")
        return int(ell)

    # This helper uses the data computed using isogeny_costs() to weigh each edge based on time complexity

    def edge_time(self, edge, p, L):
        ell = self.edge_degree(edge)
        ell_ind = L.index(ell)
        iso_time = self.path_times[str(p)][self._L_key(L)][edge[0]][ell_ind][1]
        return iso_time

    # Determine the distance of each pair in a GL graph

    def shortest_paths(self, p, L):
        """
        Determine the distance between every pair of vertices in a GL graph for a specific p using time complexity as edge weight.
        Returns:
        - The diameter of the graph,
        - List of distances of the graph.
        """

        # Form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)

        # load or create the specific area in the shortest paths data file

        self.shortest_path_times.setdefault(p_key, {})
        sp_data = self.shortest_path_times[p_key]        

        # if the data already exists, return it

        if L_key in sp_data:
            entry = sp_data[L_key]
            return entry["diam"], entry["ds"]

        # load the isogeny time complexities using the previous function

        avg_times = self.isogeny_costs(p, L)

        # edit edge_time() function so it can be used as a weight_function in .shortest_path_lengths()

        def edge_time_w(edge):
            return self.edge_time(edge, p, L)

        # use .shortest_path_lengths() with edge_time_w() as a weight function to determine the shortest paths in the GL graph

        G = get_L_graph(p, L)
        ds = []
        for v in G.vertices():
            lengths = G.shortest_path_lengths(v, by_weight = True, weight_function = edge_time_w)
            nicelengths = list(lengths.items())
            for l in nicelengths:
                ds.append(l)

        # use .diameter() to compute the diameter of the graph

        diam = G.diameter(by_weight = True, weight_function = edge_time_w)       

        # format the data

        diam_py = float(diam)
        ds_py = [[str(v), float(t)] for (v, t) in ds]

        # save and return it

        sp_data[L_key] = {"diam": diam_py, "ds": ds_py}
        self.shortest_path_times[p_key] = self._sort_L_keys(self.shortest_path_times[p_key])
        self._save_json(self.data_file2, self.shortest_path_times)

        return diam_py, ds_py

    # Now plot the data found from the above function 

    def shortest_path_time_dist(self, p, L, bin_width=0.0001, plot=True):
        """
        Compute the quantity of the distances of a GL graph within specified ranges,
        - the range, or bin width, is set to 0.0001 seconds, but can be changed,
        - if plot = True, this data will be plotted.
        Returns:
        - The distributiuon
        """

        # Form the relevant keys 

        p_key = str(int(p))
        L_key = self._L_key(L)
        L_id = "".join(str(x) for x in L)

        bin_width = float(bin_width) 

        # load or create the specific area in the analysis file

        if p_key not in self.analysis:
            self.analysis[p_key] = {"L_specific": {}, "comparisons": {}}
        else:
            self.analysis[p_key].setdefault("L_specific", {})
            self.analysis[p_key].setdefault("comparisons", {})

        Lspec = self.analysis[p_key]["L_specific"]
        Lspec.setdefault("shortest_paths_time_dist", {})

        # it is useful to make a small helper here for plotting

        def plotd(dist):
            P = None
            xs = [i * bin_width for i in range(1, len(dist) + 1)]
            pts = list(zip(xs, dist))
            P = list_plot(
                pts,
                plotjoined=True,
                marker='o',
                axes_labels=[r'$t$ (seconds)', r'count'],
            )
            filename = f"shortest_paths_time_dist_{int(p)}_{L_id}.png"
            P.save(filename)
            print(f"Saved plot to {filename}")

        # if the data already exists, return it

        if L_key in Lspec["shortest_paths_time_dist"]:
            entry = Lspec["shortest_paths_time_dist"][L_key]
            distribution = entry["distribution"]

            # and plot if requested

            if plot:
                plotd(distribution)

            return distribution

        # otherwise compute the data using .shortest_paths()

        diam, ds = self.shortest_paths(p, L)

        distances = [float(t) for (_, t) in ds]
        diam_f = float(diam)

        # sort out the number of bins based on bin_width

        nbins = int(math.ceil(diam_f / bin_width))
        nbins = max(nbins, 1)

        # and calculate the quantity of distances in each

        distribution = [0] * nbins
        for t in distances:
            i = int(t / bin_width)
            if i >= nbins:
                i = nbins - 1
            distribution[i] += 1

        # format the data 

        entry = {
            "bin_width": bin_width,
            "diam": diam_f,
            "num_pairs": int(len(distances)),
            "distribution": [int(x) for x in distribution],
        }

        # save and plot if requested

        Lspec["shortest_paths_time_dist"][L_key] = self.to_jsonable(entry)
        self._save_json(self.data_file5, self.analysis)

        if plot:
            plotd(distribution)

        return distribution


    def compare_L_dist(self, p, L1, L2, bin_width=0.0001, plot=True):
        """
        Compute shortest path distributions with specific p for multiple GL with L1 ≤ L ≤ L2 in lexicographical order and largest difference in diameter.
        - bin_width is set to 0.0001 but can be changed,
        - if plot = True, the distributions will be plotted against each other for comparison
        Returns:
        - Largest difference in diameter for all L.
        """

        # It is useful to have a helper to save the plot with the necessary legend in a consistent manner

        def save_with_forced_legend(P, filename, figsize=(12, 6)):
            fig = P.matplotlib(figsize=figsize)
            ax = fig.axes[0]
            ax.set_xlabel(r'$t$', fontsize=12)
            ax.set_ylabel(r'count', fontsize=12)
            handles, labels = ax.get_legend_handles_labels()
            if labels:
                ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
            fig.savefig(filename, bbox_inches="tight")

        # and it is useful to have a helper to build the plot

        def make_overlay_plot(chain, dist_by_L, max_diff_info, bin_width=0.001):

            color_cycle = ["red", "blue", "green", "orange", "purple",
                        "brown", "teal", "black", "magenta", "gold"]

            P = None
            max_x = 0.0
            max_y = 0

            # the plot is build by continuously adding multiple plots for each L to each other

            for idx, Ls in enumerate(chain):
                dist = dist_by_L[tuple(Ls)]
                xs = [(i + 1) * bin_width for i in range(len(dist))]
                pts = list(zip(xs, dist))

                if xs:
                    max_x = max(max_x, xs[-1])
                if dist:
                    max_y = max(max_y, max(dist))

                q = list_plot(
                    pts,
                    plotjoined=True,
                    marker='o',
                    legend_label=str(Ls),
                    color=color_cycle[idx % len(color_cycle)]
                )
                P = q if P is None else (P + q)

            xmax = max_x * 1.05 if max_x > 0 else 1.0
            ymax = max_y * 1.30 if max_y > 0 else 1.0
            P.set_axes_range(xmin=0, xmax=xmax, ymin=0, ymax=ymax)

            Li, Lj, max_d = max_diff_info
            label = f"max Δdiam = {max_d:.6g} between {Li} and {Lj}"
            P += text(label, (0.5 * xmax, 0.95 * ymax),
                    horizontal_alignment='center',
                    color='blue')

            return P

        # form the relevant keys

        p_key = str(int(p))

        # build a "chain" between L1 and L2 of Ls in lexicographical order

        A = [int(x) for x in L1]
        B = [int(x) for x in L2]

        extras = [x for x in B if x not in A]
        chain = []
        current = list(A)
        chain.append(list(current))
        for x in extras:
            current = current + [x]
            chain.append(list(current))

        # and form the key to represent this

        comp_key = f"{(self._L_key(A))}->{(self._L_key(B))}"

        # load or create the area in the analysis file

        self.analysis.setdefault(p_key, {})
        self.analysis[p_key].setdefault("L_specific", {})
        self.analysis[p_key].setdefault("comparisons", {})
        comps = self.analysis[p_key]["comparisons"]
        comps.setdefault("compare_L_dist", {})

        # initialise the dictionaries to contain the diameters and distributions 

        diam_by_L = {}
        dist_by_L = {}

        # compute the diamters and shortest path distributions for each L using .shortest_path_time_dist()

        for Ls in chain:
            dist = self.shortest_path_time_dist(p, Ls, bin_width=bin_width, plot=False)
            dist_by_L[tuple(Ls)] = dist

            diam, _ = self.shortest_paths(int(p), Ls)
            diam_by_L[tuple(Ls)] = float(diam)

        # now compute the maximum difference between diameters 

        max_diff = -1.0
        max_pair = (chain[0], chain[0])

        chain_tuples = [tuple(Ls) for Ls in chain]
        for i in range(len(chain_tuples)):
            for j in range(i + 1, len(chain_tuples)):
                Li = chain_tuples[i]
                Lj = chain_tuples[j]
                d = abs(diam_by_L[Lj] - diam_by_L[Li])
                if d > max_diff:
                    max_diff = d
                    max_pair = (list(Li), list(Lj))

        # format the data as a tuple

        max_diff_info = (max_pair[0], max_pair[1], float(max_diff))

        # format and save the data to the analysis data file

        comps["compare_L_dist"][comp_key] = {
            "L1": self._L_key(A),
            "L2": self._L_key(B),
            "max_diff_info": (self._L_key(max_pair[0]), self._L_key(max_pair[1]), float(max_diff)),
        }

        tmp = self.analysis[p_key].pop("comparisons", {})
        self.analysis[p_key]["comparisons"] = tmp

        self.analysis[p_key]["comparisons"] = dict(
            sorted(self.analysis[p_key]["comparisons"].items())
        )

        self._save_json(self.data_file5, self.analysis)

        # plot the data if requested

        def L_id(L): 
            return "".join(str(x) for x in L)

        P = None
        if plot:
            P = make_overlay_plot(chain, dist_by_L, max_diff_info, bin_width=bin_width)
            filename = f"compare_L_dist_{p}_{L_id(A)}_{L_id(B)}.png"
            save_with_forced_legend(P, filename, figsize=(12, 6))
            print(f"Saved plot to {filename}")

        return max_diff_info

    def find_recurring_dist(self, p, L1, L2):
        """
        Computes which L between L1 and L2 has shortest path distribution last to change, it can be the very last, (note this can be dependent on bin width).
        Returns:
        - Value of L,
        - Number of identical distributions.
        """

        # Form the relevant keys and make sure the data exists

        p_key = str(int(p))
        self.compare_L_dist(p, L1, L2, plot=False)

        # load the area in the analysis data file

        Lspec = self.analysis[p_key]["L_specific"]
        Lspec.setdefault("shortest_paths_time_dist", {})

        # initialise the dominant distribution as the first one until we reach a non-identical distribution

        dom_L = list(Lspec["shortest_paths_time_dist"].keys())[0]
        dom_dist = list(Lspec["shortest_paths_time_dist"].values())[0]

        # count the number of identical distributions

        identical_dist_count = 0
        for L_key in list(Lspec["shortest_paths_time_dist"].keys())[0:]:
            if Lspec["shortest_paths_time_dist"][L_key] == dom_dist:
                identical_dist_count += 1
            else:
                dom_L = L_key
                dom_dist = Lspec["shortest_paths_time_dist"][L_key]
                identical_dist_count = 0

        return dom_L, identical_dist_count

    def correlate_diam_with_eigenvals(self, p, L1, L2, Lc=None, plot=True, unramanujanness=False):
        """
        Compute the correlation between:
        - the diameter of each GL if Lc = None,
        - or the difference between the diamter of each GL and GLc,
        and 
        - the largest abs non-trivial eigenvalue of the adjacency matrix of each GL,
        - or the un-Ramanujanness of each GL (non-trivial max|λ| - 2 * sqrt(d-1)).

        - Plots the data if plot = True
        Returns:
        - The chosen values of x,
        - The chosen values of y,
        - The pearson's r correlation coefficient.
        """

        # This helper computes the pearson's r correlation coefficient

        def pearson(x, y):
            n = len(x)
            if n == 0:
                return None
            mx = sum(x) / n
            my = sum(y) / n
            vx = sum((xi - mx) ** 2 for xi in x)
            vy = sum((yi - my) ** 2 for yi in y)
            if vx == 0 or vy == 0:
                return None
            cov = sum((x[i] - mx) * (y[i] - my) for i in range(n))
            return cov / math.sqrt(vx * vy)

        # Again, it is useful to have a helper to build the plot 

        def make_scatter_plot(x_vals, y_vals, pearson_r, x_label, chain):
            color_cycle = [
                "red", "blue", "green", "orange", "purple",
                "brown", "teal", "black", "magenta", "gold"
            ]

            P = Graphics()

            for i, Ls in enumerate(chain):
                col = color_cycle[i % len(color_cycle)]
                label = self._L_key(Ls)
                P += point(
                    (x_vals[i], y_vals[i]),
                    size=35,
                    color=col,
                    legend_label=label
                )

            r_str = "undefined" if pearson_r is None else f"{pearson_r:.4f}"
            label = f"Pearson r = {r_str}"

            xmin = min(x_vals)
            xmax = max(x_vals)
            ymin = min(y_vals)
            ymax = max(y_vals)

            xrange = xmax - xmin
            yrange = ymax - ymin
            if xrange == 0:
                xrange = 1.0
            if yrange == 0:
                yrange = 1.0

            P += text(
                label,
                (xmin + 0.03 * xrange, ymax - 0.07 * yrange),
                horizontal_alignment='left',
                vertical_alignment='top',
                fontsize=11,
                color='black',
                background_color='white'
            )

            P.set_axes_range(
                xmin=xmin - 0.05 * xrange,
                xmax=xmax + 0.25 * xrange,
                ymin=ymin - 0.05 * yrange,
                ymax=ymax + 0.05 * yrange
            )

            return P

        # format L into a string for plotting purposes

        def L_id(L):
            return "".join(str(x) for x in L)

        # build the chain of Ls like before

        A = [int(x) for x in L1]
        B = [int(x) for x in L2]

        extras = [x for x in B if x not in A]
        chain = []
        current = list(A)
        chain.append(list(current))
        for x in extras:
            current = current + [x]
            chain.append(list(current))

        # form the relevant keys

        p_key = str(int(p))
        Lc_can = [int(x) for x in Lc] if (Lc is not None) else None

        comp_key = (
            f"{self._L_key(A)}->{self._L_key(B)}"
            f"|Lc={self._L_key(Lc_can) if Lc_can is not None else 'None'}"
            f"|unram={unramanujanness}"
        )

        # load or create the area in the analysis data file

        self.analysis.setdefault(p_key, {})
        self.analysis[p_key].setdefault("L_specific", {})
        self.analysis[p_key].setdefault("comparisons", {})
        comps = self.analysis[p_key]["comparisons"]
        comps.setdefault("correlate_diam_with_eigenvals", {})

        # form the relevant labels for plotting

        if unramanujanness:
            y_label = r'$\mathrm{non-trivial} $ $ \max |\lambda| - 2\sqrt{d-1}$'
        else:
            y_label = r'$\mathrm{non-trivial} $ $ \max |\lambda|$'

        x_label = r'$\mathrm{diam}$' if Lc_can is None else r'$\mathrm{diam}(L)-\mathrm{diam}(L_c)$'

        # if the data already exists return it

        if comp_key in comps["correlate_diam_with_eigenvals"]:
            entry = comps["correlate_diam_with_eigenvals"][comp_key]
            x_vals = [float(v) for v in entry["x_vals"]]
            y_vals = [float(v) for v in entry["y_vals"]]
            pearson_r = entry["pearson_r"]
            pearson_r = None if pearson_r is None else float(pearson_r)

            # and plot it if requested

            if plot:
                P = make_scatter_plot(x_vals, y_vals, pearson_r, x_label, chain)
                fname = (
                    f"correlate_diam_eigs_{p}_{L_id(A)}_{L_id(B)}_"
                    f"{('Lc_'+L_id(Lc_can)) if Lc_can is not None else 'diam'}"
                    f"{'_unram' if unramanujanness else ''}.png"
                )

                fig = P.matplotlib(figsize=(9, 6))
                ax = fig.axes[0]
                ax.set_xlabel(x_label, fontsize=12)
                ax.set_ylabel(y_label, fontsize=12)

                handles, labels = ax.get_legend_handles_labels()
                if labels:
                    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)

                fig.savefig(fname, bbox_inches="tight")
                print(f"Saved plot to {fname}")

            return x_vals, y_vals, pearson_r

        # load the eigenvalue data from the eigenvalue data file

        Lspec = self.eigenvalue_data[p_key]["L_specific"]
        Lspec.setdefault("eigenvalue_data", {})

        for Ls in chain:
            lk = self._L_key(Ls)
            if lk not in Lspec["eigenvalue_data"]:
                self.pL_eigenvalue_data(p, Ls)

        # load or compute the diameter for each L in the chain 

        diam_by_L = {}
        for Ls in chain:
            diam, _ = self.shortest_paths(p, Ls)
            diam_by_L[tuple(Ls)] = float(diam)

        # if Lc is specified, compute the diamter for Lc

        if Lc_can is not None:
            diam_c, _ = self.shortest_paths(p, Lc_can)
            diam_c = float(diam_c)

        # initialise the x and y lists

        x_vals = []
        y_vals = []

        # populate the lists with the relevant values for each L in the chain 

        for Ls in chain:
            diam_L = diam_by_L[tuple(Ls)]
            x = diam_L if Lc_can is None else (diam_L - diam_c)

            lk = self._L_key(Ls)
            nontriv_json = Lspec["eigenvalue_data"][lk]["nontrivial"]
            largest_nontriv = float(abs(nontriv_json[0]))

            if unramanujanness:
                d = sum(Ls) + len(Ls)
                threshold = float(2 * sqrt(d - 1))
                y = largest_nontriv - threshold
            else:
                y = largest_nontriv

            x_vals.append(float(x))
            y_vals.append(float(y))

        # compute the pearson's r correlation coefficient

        pearson_r = pearson(x_vals, y_vals)

        # format and save the data

        comps["correlate_diam_with_eigenvals"][comp_key] = {
            "L1": self._L_key(A),
            "L2": self._L_key(B),
            "Lc": None if Lc_can is None else self._L_key(Lc_can),
            "plot_unramanujanness": bool(unramanujanness),
            "x_vals": [float(v) for v in x_vals],
            "y_vals": [float(v) for v in y_vals],
            "pearson_r": None if pearson_r is None else float(pearson_r),
        }

        tmp = self.analysis[p_key].pop("comparisons", {})
        self.analysis[p_key]["comparisons"] = tmp
        self.analysis[p_key]["comparisons"] = dict(
            sorted(self.analysis[p_key]["comparisons"].items())
        )
        self._save_json(self.data_file5, self.analysis)

        # plot if requested

        if plot:
            P = make_scatter_plot(x_vals, y_vals, pearson_r, x_label, chain)
            fname = (
                f"correlate_diam_eigs_{p}_{L_id(A)}_{L_id(B)}_"
                f"{('Lc_'+L_id(Lc_can)) if Lc_can is not None else 'diam'}"
                f"{'_unram' if unramanujanness else ''}.png"
            )

            fig = P.matplotlib(figsize=(9, 6))
            ax = fig.axes[0]
            ax.set_xlabel(x_label, fontsize=12)
            ax.set_ylabel(y_label, fontsize=12)

            handles, labels = ax.get_legend_handles_labels()
            if labels:
                ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)

            fig.savefig(fname, bbox_inches="tight")
            print(f"Saved plot to {fname}")

        return x_vals, y_vals, pearson_r

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Experiment 3: Meet in the Middle
# -------------------------------------------------------------------------------------------------------------------------------------------------

    # This helper is useful for sorting data by increasing sensitivity

    def sort_sens_block(self, block):
        if not isinstance(block, dict):
            return block

        # average value and fixed pairs need to be handled seperately 

        avg_d_val = block.get("avg_d", None)
        fixed_pairs_data = block.get("fixed_pairs", None)
        sens_items = []
        other_items = []

        for k, v in block.items():
            if k == "avg_d":
                continue
            if k == "fixed_pairs":
                continue

            # convert the string sensitivity key into a float for size comparison

            if isinstance(k, str) and k.startswith("sensitivity_"):
                try:
                    s = float(k.split("sensitivity_", 1)[1])
                    sens_items.append((s, k, v))
                except Exception:
                    other_items.append((k, v))
            else:
                other_items.append((k, v))

        # sort the sensitivities

        sens_items.sort(key=lambda t: t[0])

        # rebuild the data block

        out = {}
        if avg_d_val is not None:
            out["avg_d"] = avg_d_val
        if fixed_pairs_data is not None:
            out["fixed_pairs"] = fixed_pairs_data
        for _, k, v in sens_items:
            out[k] = v
        for k, v in other_items:
            out[k] = v
        return out

    # Compute average pairs in terms of distance for MitM experiments

    def find_average_pair(self, p, L, sens=0.0003, recurring_index=0):
        """
        Compute the average distance in a GL graph μ, and finds the list of pairs with distance in range (μ - sens, μ + sens).
        Returns:
        - Average distance,
        - If recurring_index != None, average pair with list index recurring_index, otherwise a random choice,
        - Sensitivity used,
        - Length of the list of average pairs in this sensitivity range.
        """

        # Form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)
        sens_key = "sensitivity_" + ("{:.10g}".format(float(sens)))

        # load or create the area in the average pairs data file

        self.average_pairs.setdefault(p_key, {})
        self.average_pairs[p_key].setdefault("L_specific", {})
        self.average_pairs[p_key]["L_specific"].setdefault("average_distance_pairs", {})
        avg_pairs_block = self.average_pairs[p_key]["L_specific"]["average_distance_pairs"].setdefault(L_key, {})

        # if the data exists for this choice of sensitivity, return it

        if sens_key in avg_pairs_block and "avg_d" in avg_pairs_block:
            avg_d = avg_pairs_block["avg_d"]
            pair_list = avg_pairs_block[sens_key]
            if len(pair_list) == 0:
                print("No pairs in this sensitivity range.")
                return avg_d, None, sens, len(pair_list)

            if recurring_index is not None:
                return avg_d, tuple(pair_list[recurring_index]), sens, len(pair_list)
            return avg_d, tuple(choice(pair_list)), sens, len(pair_list)
            
        # otherwise ensure the shortest path data exists

        _, ds = self.shortest_paths(p, L)
        dist_list = [item[1] for item in ds]

        # compute the list of all the distances not equal to 0 (i.e. no (vi, vi) pairs)

        k = 1
        while True:
            if dist_list[k] == 0 or float(dist_list[k]) == 0:
                k -= 1
                break 
            k += 1
        vert_no = k
        
        dist_list_w0 = []       
        for i in range(len(dist_list)):
            if not dist_list[i] == 0 or not float(dist_list[i]) == float(0):
                dist_list_w0.append(dist_list[i])

        # use the distance list without 0 to find the mean average distance of the graph

        avg_d = sum(dist_list_w0)/len(dist_list_w0)

        # find the vertices of the graph

        verts = [item[0] for item in ds[:vert_no]]

        # initialise and populate the list of average pairs within this range 

        pair_list = []
        for i in range(vert_no):
            for j in range(vert_no):
                val = ds[i*vert_no+j][1]
                if val == 0 or float(val) == 0:
                    continue
                if val < avg_d + sens and val > avg_d - sens:
                    pair_list.append((verts[i], ds[i*vert_no+j][0], val))

        # format and save the data

        avg_pairs_block["avg_d"] = avg_d
        avg_pairs_block[sens_key] = self.to_jsonable(pair_list)

        self.average_pairs[p_key]["L_specific"]["average_distance_pairs"][L_key] = self.sort_sens_block(avg_pairs_block)

        self.average_pairs[p_key]["L_specific"]["average_distance_pairs"] = self._sort_L_keys(
            self.average_pairs[p_key]["L_specific"]["average_distance_pairs"]
        )

        self._save_json(self.data_file6, self.average_pairs)

        # return the data

        if len(pair_list) == 0:
            print("No pairs in this sensitivity range.")
            return avg_d, None, sens, len(pair_list)

        if recurring_index is not None:
            return avg_d, tuple(pair_list[recurring_index]), sens, len(pair_list)

        return avg_d, tuple(choice(pair_list)), sens, len(pair_list)

    def MitM(self, p, L, sens=0.0003, pair_index=0, fixed_pair=None):
        """
        Compute an estimate for the time complexity of a MitM attack on a GL graph between two vertices in 
        - either an average pair (in sensitivity range sens) if fixed_pair = None,
        - or a specific fixed pair, if fixed_pair != None
        Returns:
        - Total MitM complexity 
        - Shortest path time (distance)
        - Half the shortest path time (half distance)
        - Two sets of half MitM data
        """

        # This function makes use of the heapq module which is a priority queue,
        # heapq is automatically set to prioritise the smallest value in a heap,
        # perfect for prioritising neighbours of vertices for MitM

        import heapq

        # initialise the average distance and choice of pair

        avg_d = None
        pair = None

        if fixed_pair is None:
            avg_d, pair, sens, _ = self.find_average_pair(p, L, sens=sens, recurring_index=pair_index)
            if pair is None:
                return None
        else:
            pair = fixed_pair

        # initialise the vertices

        vert1 = pair[0]
        vert2 = pair[1]

        # form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)
        pair_key = f"{vert1}_&_{vert2}"
        sens_key = f"fixed_pairs"

        if fixed_pair is None:
            sens_key = "sensitivity_" + ("{:.10g}".format(float(sens)))

        # load or create the area in the MitM data file

        self.MitM_data.setdefault(p_key, {})
        self.MitM_data[p_key].setdefault("L_specific", {})
        self.MitM_data[p_key]["L_specific"].setdefault("MitM_data", {})
        MitM_block = self.MitM_data[p_key]["L_specific"]["MitM_data"].setdefault(L_key, {})

        # if the data already exists for our choice of pair, return it

        if sens_key in MitM_block:
            if pair_key in MitM_block[sens_key]:
                entry = MitM_block[sens_key][pair_key]
                return entry

        # otherwise form the graph and load or calculate the shortest path time data

        G = get_L_graph(p, L)
        ds = self.shortest_path_times[p_key][L_key]["ds"]

        # set some useful variables

        lends = len(ds)

        verts = G.vertices()
        vert_no = len(verts)

        # initialise and populate a dictionary to easily find the distances between vertices

        time_lookup = {}

        ind = 0
        knd = 0
        while True:
            row_dict = {}
            while True:
                if ds[ind][0] in row_dict.keys():
                    time_lookup[ds[knd][0]] = row_dict
                    break
                row_dict[ds[ind][0]] = ds[ind][1]
                if float(ds[ind][1]) == float(0):
                    knd = ind
                ind += 1
                if ind == lends:
                    time_lookup[ds[knd][0]] = row_dict
                    break
            if ind == lends:
                break

        # find the distance and half of the distance between the two vertices in the pair

        shortest_path = time_lookup[vert1][vert2]
        half_d = float(shortest_path) / 2.0

        # this helper reconstructs the path from a predetermined start vertex and an end vertex
        # by traversing through the parent dictionary (i.e. who is a parent of who)

        def reconstruct_path(end_vertex, parent_dict):
        
            path = []
            current = end_vertex
            while current is not None:
                path.append(current)
                current = parent_dict.get(current, None)
            path.reverse()
            return path

        # this is the half MitM function to be run twice

        def compute_half_MitM(start):
            """
            Computes half of a meet in the middle algorithm from given "start" vertex.
            Returns the relevant information: 
            - Final "winning" computed vertex
            - Final "winning" computed path and it's complexity
            - Total number of vertices claimed
            - Total complexity computed
            """

            # initialise the heap and the various dictionaries that will be useful later,
            # best_cost keeps track of the complexity of the current best path to each explored vertex,
            # best_prev_vertex keeps track of the previous vertex in the best path to each computed vertex,
            # best_len keeps track of the length of the best path to each vertex which is used for tie breaks,
            # heap is the continuously updating list which keeps track of the computed vertices and the three previously described values associated with each

            best_cost = {start: 0.0}            
            best_prev_vertex = {start: None}  
            best_len = {start: 1}    
            heap = [(0.0, start, None, 1)]
                      
            # initialise the final vertex and the complexity to it

            winner_v = None
            winner_finish = None

            # this is the main "spread"

            while heap:

                # take out the entry in the heap with the smallest total complexity,
                # to start with this will be the start vertex as just initialised 

                finish_cost, v, prev_v, path_len = heapq.heappop(heap)

                # if the complexity of the path to this vertex as recorded in this entry is greater than the complexity of the currently best recorded path, 
                # then this path is a longer path and this vertex has already been "taken", so ignore this entry

                if finish_cost != best_cost.get(v, None):
                    continue

                # since the complexity of the path is the best currently recorded complexity,
                # now check to see if the path computed has used the correct number of vertices

                if path_len != best_len.get(v, None):
                    continue

                # this is a final check for whether the entry needs to be removed

                if v != start and prev_v != best_prev_vertex.get(v, None):
                    continue

                # since the entry will always be the one with the smallest complexity, and correct path due to the checks above, if this entry has complexity greater than or equal to d/2,
                # then this entry is representative of the path between the start vertex and the middle vertex, so assign the correct values as such

                if v != start and finish_cost >= half_d:
                    winner_v = v
                    winner_finish = finish_cost
                    break

                # if the entry has complexity less than d/2, then the process not finished

                if finish_cost < half_d:

                    # for every neighbour of the vertex, find the complexity of the shortest path between them 

                    for w in G.neighbor_out_iterator(v):
                        edge_time = time_lookup[v].get(w, None)
                        if edge_time is None:
                            continue
                        edge_time = float(edge_time)
                        if w == v:
                            continue
                        if edge_time == 0.0:
                            continue

                        # for now record the total complexity to the neighbour w as the complexity to v plus the newly calculated complexity of the shortest path from v to w,
                        # it is analogous for the path length

                        new_finish = finish_cost + edge_time
                        new_len = best_len[v] + 1

                        # now check whether this path is the shortest and correct path

                        # firstly, if the vertex has not been claimed yet, claim it,
                        # secondly, if the vertex has been claimed but this path is shorter, replace it

                        eps = 1e-12
                        if (w not in best_cost) or (new_finish < best_cost[w]):
                            best_cost[w] = new_finish
                            best_prev_vertex[w] = v
                            best_len[w] = new_len
                            heapq.heappush(heap, (new_finish, w, v, new_len))

                        # thirdly, due to the nature of the shortest_paths() function, only the complexities of the shortest paths are recorded and not the paths themselves,
                        # this means the neighbour vertex w may be a neighbour through a large isogeny, 
                        # but the complexity of the path from v to w is recorded as the complexity of the shortest path which may go through other vertices, 
                        # this check ensures that if there is a previously recorded path to w with the same complexity, the path which actually traverses the relevant vertices is final

                        elif abs(new_finish - best_cost[w]) <= eps:
                            if new_len > best_len[w]:
                                best_prev_vertex[w] = v
                                best_len[w] = new_len
                                heapq.heappush(heap, (new_finish, w, v, new_len))

            # use reconstruct_path() to compute the winning path

            winning_path = None
            if winner_v is not None:
                winning_path = reconstruct_path(winner_v, best_prev_vertex)

            # now determine which paths should contribute to the total complexity

            # firstly, determine which vertices claimed spread to other vertices

            final_spread_vertices = set()
            for child, parent in best_prev_vertex.items():
                if parent is not None:
                    final_spread_vertices.add(parent)

            # initialise a set of endpoint vertices which may contribute and add the winner vertex

            contributor_endpoints = set()

            if winner_v is not None:
                contributor_endpoints.add(winner_v)

            # secondly, determine which vertices are the parents of vertices with path time complexity greater than d/2, 
            # and which are dead end vertices

            for u, c in best_cost.items():
                if u == winner_v:
                    continue

                if c >= half_d:
                    pv = best_prev_vertex.get(u, None)
                    if pv is not None:
                        contributor_endpoints.add(pv)

                elif u not in final_spread_vertices:
                    contributor_endpoints.add(u)

            # to avoid double counting, keep track of the edges and vertices counted

            used_edges = set()
            used_vertices = set()

            # now, starting at each end point, recreate the path back to the starting vertex, 
            # adding each computed edge only once,
            # this avoids double counting incidents like two > d/2 vertices having the same parent vertex,
            # or such a parent vertex also leading to a dead, among others

            for endpoint in contributor_endpoints:
                curr = endpoint
                while curr is not None:
                    parent = best_prev_vertex.get(curr, None)
                    used_vertices.add(curr)

                    if parent is not None:
                        used_edges.add((parent, curr))

                    curr = parent

            # add the complexity of each individual edge to the total complexity once

            total_comp = 0.0
            for parent, child in used_edges:
                total_comp += float(best_cost[child] - best_cost[parent])

            # format the information

            info = {
                "winner_vertex": str(winner_v),
                "winner_finish": float(winner_finish),
                "winning_path": [str(v) for v in winning_path],
                "no_vertices_claimed": int(len(used_vertices)),
                "total_comp": float(total_comp)
            }
            return float(total_comp), info

        # compute half MitM from both vertices

        total1, info1 = compute_half_MitM(vert1)
        total2, info2 = compute_half_MitM(vert2)

        # format and save the information

        entry = {
            "total_comp": total1+total2,
            "shortest_path_time": shortest_path,
            "half_shortest_path_time": half_d,
            str(vert1): info1,
            str(vert2): info2,
        }

        MitM_block.setdefault(sens_key, {})
        MitM_block[sens_key][pair_key] = self.to_jsonable(entry)

        self.MitM_data[p_key]["L_specific"]["MitM_data"][L_key] = self.sort_sens_block(MitM_block)

        self.MitM_data[p_key]["L_specific"]["MitM_data"] = self._sort_L_keys(
            self.MitM_data[p_key]["L_specific"]["MitM_data"]
        )

        self._save_json(self.data_file4, self.MitM_data)

        return entry

    # Now statistical analysis of MitM

    def compute_MitM_spread_for_pL(self, p, L, fixed_sens=None, fixed_pairs=None):
        """
        Compute the mean and standard deviation of the MitM complexities of a number of average pairs in a GL graph,
        - if fixed_sens = None, then the sensitivity is the smallest such that contains at least 8 pairs, 
        otherwise fixed_sens must contain at least 8 pairs,
        - or specify the exact pairs (thus any number of pairs) using fixed_pairs.
        Returns:
        - Average MitM complexity,
        - Standard deviation,
        - Min MitM complexity,
        - Max MitM complexity,
        - Number of pairs used,
        - Sensitivity. 
        """

        # Form the relevant keys

        p_key = str(int(p))
        L_key = self._L_key(L)
        fixed_sens_key = "sensitivity_" + ("{:.10g}".format(float(fixed_sens))) if fixed_sens is not None else None

        # load or create the area in the analysis data file

        self.analysis.setdefault(p_key, {})
        self.analysis[p_key].setdefault("L_specific", {})
        self.analysis[p_key]["L_specific"].setdefault("MitM_data", {})
        MitM_comp_data = self.analysis[p_key]["L_specific"]["MitM_data"].setdefault("stats", {})

        # if the sensitivity has been specified, check the data already exists

        if fixed_sens is not None:
            if fixed_sens_key in MitM_comp_data:
                if L_key in MitM_comp_data[fixed_sens_key]:
                    entry = MitM_comp_data[fixed_sens_key][L_key]
                    return entry

        # define the minimum number of pairs, 8 is arbitrary, this can be changed

        min_pair_no = 8

        # initialise the list for all the MitM data and the pair number variable

        total_comp_list = []
        pair_no = 0

        sens = "fixed_pairs"
        sens_key = "fixed_pairs"

        # if the pairs are controlled, the function is different

        if fixed_pairs is None:      

            # determine the right sensitivity range to include min_pair_no number of pairs

            # increments of 0.0000001 has proven to be small enough to gradually encapsulate an appropriate amount of pairs,
            # instead of large numbers all at once due to the peaks seen in the shortest path distributions

            sens = 0.0000001

            if fixed_sens == None:
                while pair_no < min_pair_no:
                    _1, _2, _3, pair_no = self.find_average_pair(p, L, sens=sens)
                    sens += 0.0000001

                sens = sens-0.0000001

            # if fixed_sens doesn't encapsulate enough pairs, the function will say so

            else:
                sens = fixed_sens
                _1, _2, _3, pair_no = self.find_average_pair(p, L, sens=sens)
                if pair_no < min_pair_no:
                    print(f"Less than {min_pair_no} pairs in this sensitivity range.")
                    return None

            sens_key = "sensitivity_" + ("{:.10g}".format(float(sens)))

            # check if the data already exists in this sensitivity range

            if sens_key in MitM_comp_data:
                if L_key in MitM_comp_data[sens_key]:
                    entry = MitM_comp_data[sens_key][L_key]
                    return entry

            # start populating the list of data

            for i in range(pair_no):
                info = self.MitM(p, L, sens=sens, pair_index=i)
                total_comp_list.append(info["total_comp"])
        
        # if the pairs are predetermined, this function is a lot less complicated

        else:
            pair_no = len(fixed_pairs)
            for i in range(pair_no):
                info = self.MitM(p, L, fixed_pair=fixed_pairs[i])
                total_comp_list.append(info["total_comp"])

        # compute the relevant data

        avg_comp = statistics.mean(total_comp_list)
        stdv = statistics.stdev(total_comp_list)
        max_val, min_val = max(total_comp_list), min(total_comp_list)

        # format, save, and return it

        entry = {
            "average_complexity": avg_comp,
            "standard_deviation": stdv,
            "maximum": max_val,
            "minimum": min_val,
            "pair_number": pair_no,
            "sensitivity": sens
        }

        MitM_comp_data.setdefault(sens_key, {})
        MitM_comp_data[sens_key][L_key] = self.to_jsonable(entry)

        self.analysis[p_key]["L_specific"]["MitM_data"]["stats"][sens_key] = self._sort_L_keys(
            self.analysis[p_key]["L_specific"]["MitM_data"]["stats"][sens_key]
        )

        self.analysis[p_key]["L_specific"]["MitM_data"]["stats"] = self.sort_sens_block(MitM_comp_data)

        self._save_json(self.data_file5, self.analysis)

        return entry
            
    # Now plot and visualise the MitM data

    def compute_and_plot_MitM_data(self, p, L1, L2, fixed_sens=None, fix_pairs=False, fixed_no_pairs=8, plot=True):
        """
        Compute MitM satistical data and un-Ramanujanness with specific p for multiple GL with L1 ≤ L ≤ L2 in lexicographical order,      
        - for the average pairs for each graph (8 minimum) unless fix_pairs = True, in which case the pairs are fixed to be the ones 
        from the GL graph with the minimum number of average pairs, which has a fixed minimum of fixed_no_pairs, 
        - if fix_pairs = False, optionally change the fixed_sens,
        - if plot = True, MitM data is plotted against un-Ramanujanness.
        Returns (for each GL graph):
        - Un-Ramanujanness,
        - Average MitM,
        - Standard deviation,
        - Number of pairs,
        and
        - Dominant shortest path distribution (dominant L),
        - Sensitivity.
        """

        # Form the relevant keys

        p_key = str(int(p))
        chain_L_key = f"{(self._L_key(L1))}->{(self._L_key(L2))}"

        # define the minimum pair number unless pairs are fixed

        min_pair_no = 8
        if fix_pairs:
            chain_L_key += chain_L_key + " (fixed pairs)"
            min_pair_no = fixed_no_pairs

        # load or create the area in the analysis data file

        self.analysis.setdefault(p_key, {})
        self.analysis[p_key].setdefault("comparisons", {})
        MitM_comparisons_block = self.analysis[p_key]["comparisons"].setdefault("MitM_comparisons", {})

        # define the relevant lists of data to be populated

        eigenval_diffs = []
        average_complexities = []
        standard_deviations = []
        pair_numbers = []
        dom_L = None
        sens = None

        # form the lexicographical chain of Ls 

        A = [int(x) for x in L1]
        B = [int(x) for x in L2]

        extras = [x for x in B if x not in A]
        chain = []
        current = list(A)
        chain.append(list(current))
        for x in extras:
            current = current + [x]
            chain.append(list(current))

        # if the data already exists, then no computation is required

        if chain_L_key in MitM_comparisons_block:
            entry = MitM_comparisons_block[chain_L_key]
            for Ls in list(entry.keys())[:-2]:
                eigenval_diffs.append(entry[Ls]["eigenval_threshold_diff"])
                average_complexities.append(entry[Ls]["average_MitM"])
                standard_deviations.append(entry[Ls]["standard_deviation"])
                pair_numbers.append(entry[Ls]["pair_number"])
            dom_L = entry["dominant_dist"]
            sens = entry["sensitivity"]

        # now determine the sensitivity to encapsulate the minimum number of pairs

        else:
            min_sens = 0.0000001
            pair_no_list = []
            for Ls in chain:
                _1, _2, _3, pair_no = self.find_average_pair(p, Ls, sens=min_sens)
                while pair_no < min_pair_no:
                    min_sens += 0.0000001
                    _1, _2, _3, pair_no = self.find_average_pair(p, Ls, sens=min_sens)
                pair_no_list.append(pair_no)

            # make note of which L gives the minimum number of pairs in case fix_pairs = True 

            min_pair_no = min(pair_no_list)
            chain_ind = pair_no_list.index(min_pair_no)
            min_L = chain[chain_ind]

            # if fixed_sens is specified, make sure fixed_sens is greater than or equal to min_sens

            sens = min_sens
            if fixed_sens is not None:
                if fixed_sens > min_sens:
                    sens = fixed_sens
                elif fixed_sens < min_sens:
                    print("Sensitivity is too small; increase or leave blank.")
                    return None

            # if fix_pairs = True, define the necessary pairs 

            fixed_pairs = []
            if fix_pairs:
                for i in range(min_pair_no):
                    _1, pair, _3, _4 = self.find_average_pair(p, min_L, sens=sens, recurring_index=i)
                    fixed_pairs.append(pair)

            # now for each L, populate the data lists 

            for Ls in chain:

                # with un-Ramanujanness

                eigenval_data = self.pL_eigenvalue_data(p, Ls)
                largest_eigenval = eigenval_data["nontrivial"][0]
                d = sum(Ls) + len(Ls)
                threshold = 2 * sqrt(d - 1)
                diff = abs(largest_eigenval) - threshold
                eigenval_diffs.append(diff)

                # and MitM data using .compute_MitM_spread_for_pL()

                Ls_entry = None
                if fix_pairs:
                    Ls_entry = self.compute_MitM_spread_for_pL(p, Ls, fixed_pairs=fixed_pairs)
                else:
                    Ls_entry = self.compute_MitM_spread_for_pL(p, Ls, fixed_sens=sens)
                avg_comp = Ls_entry["average_complexity"]
                stdv = Ls_entry["standard_deviation"]
                pair_number = Ls_entry["pair_number"]
                average_complexities.append(avg_comp)
                standard_deviations.append(stdv)
                pair_numbers.append(pair_number)

            # copmute the dominant distribution for the sake of analysis

            dom_L, _ = self.find_recurring_dist(p, L1, L2)

        # plot if requested

        if plot:
            color_cycle = [
                "red", "blue", "green", "orange", "purple",
                "brown", "teal", "black", "magenta", "gold"
            ]
            
            def fmt_sens(x):
                return f"{float(x):.6g}"

            P = Graphics()

            xs = [float(v) for v in eigenval_diffs]
            ys = [float(v) for v in average_complexities]
            ss = [float(v) for v in standard_deviations]

            xmin, xmax = min(xs), max(xs)
            ymin = min(y - s for (y, s) in zip(ys, ss))
            ymax = max(y + s for (y, s) in zip(ys, ss))

            xrange = xmax - xmin
            yrange = ymax - ymin
            if xrange == 0: xrange = 1.0
            if yrange == 0: yrange = 1.0

            pad_x = 0.05 * xrange
            pad_y = 0.08 * yrange

            cap = 0.01 * xrange

            for i, Ls in enumerate(chain):
                x = float(eigenval_diffs[i])
                y = float(average_complexities[i])
                s = float(standard_deviations[i])
                col = color_cycle[i % len(color_cycle)]

                q = point((x, y), size=35, color=col, legend_label=self._L_key(Ls) + f"  {pair_numbers[i]} pairs")

                ebar = line([(x, y - s), (x, y + s)], color=col, thickness=1.2)
                caps = line([(x - cap, y + s), (x + cap, y + s)], color=col, thickness=1.2) + \
                    line([(x - cap, y - s), (x + cap, y - s)], color=col, thickness=1.2)

                P += (q + ebar + caps)

            x_dummy = xmax + 10 * pad_x
            y_dummy = ymax + 10 * pad_y

            P += point((x_dummy, y_dummy), size=0, color="white",
                    legend_label=f"dominant dist = {dom_L}")
            P += point((x_dummy, y_dummy), size=0, color="white",
                    legend_label=f"sensitivity = {fmt_sens(sens)}")
            if fix_pairs:
                P += point((x_dummy, y_dummy), size=0, color="white",
                        legend_label=f"fixed pairs")

            P.set_axes_range(
                xmin=xmin - pad_x,
                xmax=xmax + pad_x + 0.30 * xrange,
                ymin=ymin - pad_y,
                ymax=ymax + pad_y
            )

            filename = f"MitM_scatter_{int(p)}_{''.join(map(str, A))}_{''.join(map(str, B))}.png"
            fig = P.matplotlib(figsize=(10, 6))
            ax = fig.axes[0]
            handles, labels = ax.get_legend_handles_labels()
            if labels:
                ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
            ax.set_xlabel(r'$\mathrm{non-trivial} $ $ \max |\lambda| - 2\sqrt{d-1}$', fontsize=12)
            ax.set_ylabel("Avg MitM time complexity", fontsize=12)

            fig.savefig(filename, bbox_inches="tight")
            print(f"Saved plot to {filename}")

        # return the data if it already exists

        if chain_L_key in MitM_comparisons_block:
            return MitM_comparisons_block[chain_L_key]

        # format, save, and return the data

        entry = {

            }

        for i, Ls in enumerate(chain):
            entry[self._L_key(Ls)] = {
                "eigenval_threshold_diff": eigenval_diffs[i],
                "average_MitM": average_complexities[i],
                "standard_deviation": standard_deviations[i],
                "pair_number": pair_numbers[i]
                }

        entry["dominant_dist"] = dom_L
        entry["sensitivity"] = sens

        MitM_comparisons_block[chain_L_key] = entry
        self._save_json(self.data_file5, self.analysis)

        return entry

# -------------------------------------------------------------------------------------------------------------------------------------------------
# References
# -------------------------------------------------------------------------------------------------------------------------------------------------

# [1] https://eprint.iacr.org/2025/155.pdf
# [2] https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/isogeny_small_degree.html#sage.schemes.elliptic_curves.isogeny_small_degree.isogenies_prime_degree

# =================================================================================================================================================

