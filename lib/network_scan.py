def likelyhood_ratio(baseline, observed):
    if observed < 0:
        raise ValueError("The number of observed events is negative")
    if baseline < 0:
        raise ValueError("The number of simulated events is negative")
    if baseline == 0 and observed > 0:
        return 1.0
    if observed > baseline:
        return np.power((observed / baseline), observed) * np.exp(baseline - observed)
    return 1.0
def continuous_homogeneous_possion_likelyhood_test(baseline, observed, total_baseline, total_observed):
    if observed < 0:
        raise ValueError("The number of observed events is negative")
    if baseline < 0:
        raise ValueError("The number of simulated events is negative")
    if observed > baseline:
        return (np.power((observed / baseline), observed) * np.power((total_observed - observed)/ (total_baseline - baseline), total_observed - observed)) / np.power((total_observed/total_baseline), total_observed)
    
    return 1.0
    
    
def expected_counts(study_G_size, G_size, G_counts):
    return G_counts * study_G_size / G_size
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
            count += int(G.nodes[p]['event_count'])
    return count

def network_scan_statistic_at_point(dic, pi, G, wholeN, L, targetSize, step = 20, start = 10000, max_distance = 25000, global_HS = True, plot = False, show_log = False):
    if(show_log):
        print(pi)
    time1 = time.time()
    if(show_log):
        print("Compute extended-tree")
    et = extended_shortest_path_tree_dict(pi, G, cutoff = max_distance) 
    LT = Ltp(et, max_distance)
    while(LT < targetSize):
        start = max_distance
        max_distance *= 2
        et = extended_shortest_path_tree_dict(pi, G, cutoff = max_distance) 
        LT_double = Ltp(et, max_distance)
        if not LT_double > LT:
            dic[pi] = 1.0
            if(show_log):
                print(pi, " on island with size ", LT_double,  ", smaller than targetSize ", targetSize)
            return
        LT = LT_double
        if(show_log):
            print("Double the etree max_distance = ", max_distance, " new tree size: ", LT_double)
        
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
    
    t = start
    l = Ltp(et, t)
    while(l < targetSize):
        t += step
        l = Ltp(et, t)
    if(show_log):
        print("Subnet size: ", l)

    expected_n = wholeN / L * l 
    observed_n = ntpi(t, pi, G,matched_index)
    likelyhood = continuous_homogeneous_possion_likelyhood_test(expected_n, observed_n, wholeN, wholeN)
    if(show_log):
        print("Expected: ", expected_n, " Observed: ", observed_n, " likelyhood at ", pi, " is ", likelyhood)
    #return likelyhood
    dic[pi] = likelyhood
    return likelyhood

def network_scan_on_G(G, center_points, no_of_processes, window_size = 0):
    p = multiprocessing.Pool(no_of_processes)
    d = multiprocessing.Manager().dict()
    process_list = []
    N = wholeN(G)
    LG = L(G)
    if(window_size == 0):
        targetSize = LG/10
    else:
        target_size = window_size
    for pi in center_points:
        process_list.append((d, pi, G, N, LG, targetSize))
    p.starmap(network_scan_statistic_at_point, process_list)
    return d