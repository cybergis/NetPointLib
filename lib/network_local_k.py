from hotspot import *

def Var(G, t, l, pi, n, et):
    Ltpi = Ltp(et, t)
    if n == 0:
        n = 1
    var = 1/ n * l * Ltpi * (1-Ltpi/l)
    return var

def Ktp_grouped(G, t, l, pi, n, matched_index):
    if n == 0:
        n = 1
    return l/(n)*ntpi_grouped(t, pi, G, matched_index)

def ntpi_grouped(t, pi, G,matched_index):
    count = 0
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=t, weight = 'length') #No need to start from scratch every time
    for p in paths.keys():
        if p in matched_index and (p not in found):
            #print(p)
            #print(p in found)
            found.append(p)
            if 'event_count' in G.nodes[p]:
                count += int(G.nodes[p]['event_count'])
            else:
                count += 1
    return count

def network_local_k(pi, G, wholeN, l, alpha = 0.005, t = 1600, m = 1, global_HS = True, plot = False, show_log = True):
    if(show_log):
        print(pi)
    time1 = time.time()
    if(show_log):
        print("Compute extended-tree")
    et = extended_shortest_path_tree_dict(pi, G, cutoff = t*m) 
    matched_index = list()
    if(show_log):
        print("Start dijkstra's")
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=t, weight = 'length') #No need to start from scratch every time
    for p in paths.keys():
        if p < 0:
            matched_index.append(p)
    
    #pi = matched_index[center_index]
    cluster = 0
    if(show_log):
        timePrep = time.time()
        print("Time for preparation: ", str(timePrep - time1))
    file1 = open("cluster.txt", "a")  # append mode
    #file1.write(str(prob)+ ", ")
    if not global_HS:
        H = G.subgraph(list(paths.keys()))
        G = H
    for i in range(1, m+1):
        if show_log:
            print("i = ", i)
            print("No. of points in range: ", ntpi(t*i, pi, G ,matched_index))
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
        std = Var(G, t * i, l, pi, wholeN, et)**0.5
        sTtime = time.time()
        if show_log:
            print("std cal: ", std)
        k = Ktp_grouped(G, t * i, l, pi, wholeN, matched_index)
        ktime = time.time()
        if show_log:
            print("k cal: ",k)
        normalDist = scipy.stats.norm(mean, std)
        if k >= normalDist.ppf(1-alpha):
            if show_log:
                print("At ti = ", t*i, " k >upper critical")
            cluster = cluster + 1
            #dic[n] = 1
            return 1
        elif k < normalDist.ppf(alpha):
            if show_log:
                print("At ti = ", t*i, " k <lower critical")
            #dic[n] = -1
            return -1
        else:
            if show_log:
                print("At ti = ", t*i, " CSR")
            #dic[n] = 0
            return 0
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
        #print("Saving result for " + str(prob) + "to " + "./" + classPath + "/" + str(c_code) + "_" + str(t) + "_" + str(prob))
    if show_log:
        print("Total count: ", str(cluster))
    file1.write(str(cluster) + ", "+ str(time2 - time1) + "\n")
    file1.close()
    #remove_node(G, center_index)
    return G