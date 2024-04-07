def detectHotspot_matched_prob(prob, G, t=20, m=80, c_code="default", max_distance=4800, global_HS=False, plot=True, show_log=False):
    """
    Detects hotspots in a specified area based on the graph G and probability thresholds.

    This function identifies hotspots by analyzing the distribution of matched events within the graph G, considering
    the spatial probability and temporal thresholds. It can operate in a global or localized context, modifying
    the graph as necessary to highlight hotspots and potentially de-emphasize non-hotspots areas.

    Args:
        prob (tuple): The longitude and latitude coordinates representing the center of analysis.
        G (networkx.MultiDiGraph): The directed graph representing the street network for hotspot analysis.
        t (int, optional): Temporal threshold for hotspot analysis. Defaults to 20.
        m (int, optional): Spatial multiplier for expanding the analysis area. Defaults to 80.
        c_code (str, optional): Category code or identifier for the analysis. Defaults to "default".
        max_distance (int, optional): Maximum distance for considering points in the hotspot analysis. Defaults to 4800 meters.
        global_HS (bool, optional): Flag to determine if the hotspot analysis should be global. Defaults to False.
        plot (bool, optional): Flag to enable or disable plotting of the hotspot analysis results. Defaults to True.
        show_log (bool, optional): Flag to enable or disable logging of the process details. Defaults to False.

    Returns:
        networkx.MultiDiGraph: The modified graph G, with hotspots highlighted based on the analysis.
    """
    if(show_log):
        print(prob)
    
    pi = ox.distance.nearest_nodes(G, prob[0], prob[1])
    if(show_log):
        print(pi)
    time1 = time.time()
    if(show_log):
        print("Compute extended-tree")
    et = extended_shortest_path_tree(pi, G, cutoff = t*m) 
    matched_index = list()
    if(show_log):
        print("Start dijkstra's")
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=max_distance, weight = 'length') #No need to start from scratch every time
    for p in paths.keys():
        if p < 0:
            matched_index.append(p)
    
    #pi = matched_index[center_index]
    cluster = 0
    if(show_log):
        timePrep = time.time()
        print("Time for preparation: ", str(timePrep - time1))
    file1 = open("cluster.txt", "a")  # append mode
    file1.write(str(prob)+ ", ")
    if not global_HS:
        H = G.subgraph(list(paths.keys()))
        G = H
    for i in range(1, m+1):
        if show_log:
            print("i = ", i)
            print("No. of points in range: ", ntpi(t*i, pi, G,matched_index))
        if ntpi(t*i, pi, G,matched_index) == 0:
            if show_log:
                print("No event in range, skip")
            cluster = cluster - 1
            continue
        stime = time.time()
        mean = Ltp(et, t * i)
        ltime = time.time()
        if show_log:
            print("Ltp cal: ", mean)
        std = Var(G, t * i, pi, wholeN(G), et)**0.5
        sTtime = time.time()
        if show_log:
            print("std cal: ", std)
        k = Ktp(G, t * i, pi, wholeN(G), matched_index)
        ktime = time.time()
        if show_log:
            print("k cal: ",k)
        normalDist = scipy.stats.norm(mean, std)
        if k >= normalDist.ppf(1-.01):
            file1.write("1, ")
            if show_log:
                print("At ti = ", t*i, " k >upper critical")
            cluster = cluster + 1
        elif k < normalDist.ppf(.01):
            file1.write("-1, ")
            if show_log:
                print("At ti = ", t*i, " k <lower critical")
            cluster = cluster - 1
        else:
            file1.write("0, ")
            if show_log:
                print("At ti = ", t*i, " CSR")
            #cluster = cluster - 1
    time2 = time.time()
    if plot:
        nc = ["r" if (node <0) else "b" for node in G.nodes()]
        ns = [10 if (node <0) else 0 for node in G.nodes()]
        classPath = str(c_code)+"_" +str(t) + "_" + str(m)+"/"
        if cluster > 0:
            classPath  = classPath + "Cluster"
        elif cluster < 0:
            classPath = classPath + "De_cluster"
        else:
            classPath = classPath + "CSR"
        fig, ax = ox.plot_graph(G, node_color=nc, node_size = ns, bgcolor = '#ffffff', save = True, filepath = "./" + classPath + "/" + str(c_code) + "_" + str(t) + "_" + str(prob) + ".jpg")
        print("Saving result for " + str(prob) + "to " + "./" + classPath + "/" + str(c_code) + "_" + str(t) + "_" + str(prob))
    if show_log:
        print("Total count: ", str(cluster))
    file1.write(str(cluster) + ", "+ str(time2 - time1) + "\n")
    file1.close()
    #remove_node(G, center_index)
    return G