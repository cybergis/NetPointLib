import osmnx as ox
import random
import networkx as nx
import numpy as np
import pandas as pd
import time
import math
import shapely
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform, pdist
from sklearn.cluster import DBSCAN
import warnings
import scipy
import geopy.distance
from anytree import Node, RenderTree, AsciiStyle, find, PreOrderIter
warnings.filterwarnings('ignore')
ox.config(use_cache=True, log_console=False)
def pedal(p1, p2, p3):
    if p2[0] != p1[0]:
        #print(p1, p2, p3)
        k, b = np.linalg.solve([[p1[0], 1], [p2[0], 1]], [p1[1], p2[1]])
        x = np.divide(((p2[0] - p1[0]) * p3[0] + (p2[1] - p1[1]) * p3[1] - b * (p2[1] - p1[1])),
                      (p2[0] - p1[0] + k * (p2[1] - p1[1])))
        y = k * x + b

    else:
        x = p1[0]
        y = p3[1]
    d1 = (x - p1[0])**2 + (y - p1[1])**2
    d2 = (x - p2[0])**2 + (y - p2[1])**2
    d = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2
    if d1 > d:
        return p2[0], p2[1]
    elif d2 > d:
        return p1[0], p1[1]
    else:
        return x, y
def fix_oneway(G):
    #edge_oneway = nx.get_edge_attributes(G, "oneway")
    for e in G.edges:
        if G.edges[e]['oneway']:
            G.add_edge(e[1], e[0], 0)
            attrs = {(e[1], e[0], 0):{"length": G.edges[e]['length'], "oneway": False}, e:{"length": G.edges[e]['length'], "oneway": False}}
            nx.set_edge_attributes(G, attrs)
def match_event_to_edges_grouped(events, G, group_dist = 10, show_mod = False): #May need to create a function for reverting the ops
    #max_index = max(G.nodes) + 1
    #max_index = base_id + 1
    matched_points = []
    for i in range(len(events)):
        #event_index = -1 * (max_index + i)
        event_index = events.iloc[i].name * -1
        nearest_node = ox.nearest_nodes(G, events.iloc[i]['x'], events.iloc[i]['y'])
        
        distance_to_nearest_node = geopy.distance.geodesic((events.iloc[i]['y'], events.iloc[i]['x']), (G.nodes[nearest_node]['y'], G.nodes[nearest_node]['x'])).m
        if(show_mod):
            print("Nearest node " + str(nearest_node) + " " + str(distance_to_nearest_node))
        if (nearest_node < 0 and distance_to_nearest_node < group_dist):
            G.nodes[nearest_node]['event_count'] += 1
            if(show_mod):
                print("Event count of node " + str(nearest_node) + " has increased to " + str(G.nodes[nearest_node]['event_count']))
        else:
            if (show_mod):
                print("Finding nearest edge for ", event_index, ", ", events.iloc[i])
            nearest_edge = ox.nearest_edges(G, events.iloc[i]['x'], events.iloc[i]['y'])
            if (show_mod):
                print("Nearest edge is ", nearest_edge)
            p1 = nearest_edge[0]
            p2 = nearest_edge[1]
            if(euclidean_distance(G,p1, p2) < 1):
                if(p1 < 0):
                    G.nodes[p1]['event_count'] += 1
                    if(show_mod):
                        print("Short edge: Event count of node " + str(p1) + " has increased to " + str(G.nodes[nearest_node]['event_count']))
                    return
                elif(p2 < 0):
                    G.nodes[p2]['event_count'] += 1
                    if(show_mod):
                        print("Short edge: Event count of node " + str(p2) + " has increased to " + str(G.nodes[nearest_node]['event_count']))
                    return
            try:
                match_x, match_y = pedal((G.nodes[p1]['x'], G.nodes[p1]['y']), (G.nodes[p2]['x'], G.nodes[p2]['y']), (events.iloc[i]['x'], events.iloc[i]['y']))
            except Exception as e:
                raise ValueError('Exception ' + str(e) + ' wheh computing' +str((G.nodes[p1]['x'], G.nodes[p1]['y'])) + str ((G.nodes[p2]['x'], G.nodes[p2]['y'])) + str ((events.iloc[i]['x'], events.iloc[i]['y'])))
            if(show_mod):
                print("Match ", match_x, match_y)
            
            #event_index = -1 * (max_index + i)
            edge_lengths = nx.get_edge_attributes(G, "length")
            originalLength = edge_lengths[(p1, p2, 0)]
            x_dif = G.nodes[p2]['x'] - G.nodes[p1]['x']
            y_dif = G.nodes[p2]['y'] - G.nodes[p1]['y']
            if not y_dif == 0:
                length_1e = (match_y - G.nodes[p1]['y'])/ (G.nodes[p2]['y'] - G.nodes[p1]['y']) * originalLength
            elif not x_dif== 0:
                length_1e = (match_x - G.nodes[p1]['x'])/ (G.nodes[p2]['x'] - G.nodes[p1]['x']) * originalLength
            else:
                length_1e = 0
            if length_1e < 1e-7:
                length_1e = 0
            length_2e = originalLength - length_1e
            if length_2e < 1e-7:
                length_2e = 0
            G.add_node(event_index)
            G.nodes[event_index]['x'] = match_x
            G.nodes[event_index]['y'] = match_y
            G.nodes[event_index]['street_count'] = 2
            G.nodes[event_index]['event_count'] = 1
            G.remove_edges_from(((p1, p2, 0), (p2, p1, 0)))
            if(show_mod):
                print("Remove ", p1, "-" , p2)
                print("Add ", p1, "-", event_index, "-", p2)
            G.add_edges_from(((p1, event_index, 0), (event_index, p1, 0), (p2, event_index, 0), (event_index, p2, 0)))
            attrs = {(p1, event_index, 0):{"length": length_1e, "oneway": False}, (event_index, p1, 0):{"length": length_1e, "oneway": False}, (p2, event_index, 0):{"length": length_2e, "oneway": False}, (event_index, p2, 0):{"length": length_2e, "oneway": False}}
            nx.set_edge_attributes(G, attrs)
        matched_points.append(event_index)
    return matched_points
def match_event_to_edges(events, G, show_mod = False): #May need to create a function for reverting the ops
    #max_index = max(G.nodes) + 1
    #max_index = base_id + 1
    matched_points = []
    for i in range(len(events)):
        #event_index = -1 * (max_index + i)
        event_index = events.iloc[i].name * -1
        if (show_mod):
            print("Finding nearest edge for ", event_index, ", ", events.iloc[i])
        nearest_edge = ox.nearest_edges(G, events.iloc[i]['x'], events.iloc[i]['y'])
        if (show_mod):
            print("Nearest edge is ", nearest_edge)
        p1 = nearest_edge[0]
        p2 = nearest_edge[1]
        match_x, match_y = pedal((G.nodes[p1]['x'], G.nodes[p1]['y']), (G.nodes[p2]['x'], G.nodes[p2]['y']), (events.iloc[i]['x'], events.iloc[i]['y']))
        if(show_mod):
            print("Match ", match_x, match_y)
        
        #event_index = -1 * (max_index + i)
        edge_lengths = nx.get_edge_attributes(G, "length")
        originalLength = edge_lengths[(p1, p2, 0)]
        x_dif = G.nodes[p2]['x'] - G.nodes[p1]['x']
        y_dif = G.nodes[p2]['y'] - G.nodes[p1]['y']
        if not y_dif == 0:
            length_1e = (match_y - G.nodes[p1]['y'])/ (G.nodes[p2]['y'] - G.nodes[p1]['y']) * originalLength
        elif not x_dif== 0:
            length_1e = (match_x - G.nodes[p1]['x'])/ (G.nodes[p2]['x'] - G.nodes[p1]['x']) * originalLength
        else:
            length_1e = 0
        if length_1e < 1e-7:
            length_1e = 0
        length_2e = originalLength - length_1e
        if length_2e < 1e-7:
            length_2e = 0
        G.add_node(event_index)
        G.nodes[event_index]['x'] = match_x
        G.nodes[event_index]['y'] = match_y
        G.nodes[event_index]['street_count'] = 2
        G.remove_edges_from(((p1, p2, 0), (p2, p1, 0)))
        if(show_mod):
            print("Remove ", p1, "-" , p2)
            print("Add ", p1, "-", event_index, "-", p2)
        G.add_edges_from(((p1, event_index, 0), (event_index, p1, 0), (p2, event_index, 0), (event_index, p2, 0)))
        attrs = {(p1, event_index, 0):{"length": length_1e, "oneway": False}, (event_index, p1, 0):{"length": length_1e, "oneway": False}, (p2, event_index, 0):{"length": length_2e, "oneway": False}, (event_index, p2, 0):{"length": length_2e, "oneway": False}}
        nx.set_edge_attributes(G, attrs)
        matched_points.append(event_index)
    return matched_points

def match_event_to_edges_distance(events, G, show_mod = False): #May need to create a function for reverting the ops
    #max_index = max(G.nodes) + 1
    #max_index = base_id + 1
    matched_points = []
    for i in range(len(events)):
        #event_index = -1 * (max_index + i)
        event_index = events.iloc[i].name * -1
        if (show_mod):
            print("Finding nearest edge for ", event_index, ", ", events.iloc[i])
        nearest_edge = ox.nearest_edges(G, events.iloc[i]['x'], events.iloc[i]['y'])
        if (show_mod):
            print("Nearest edge is ", nearest_edge)
        p1 = nearest_edge[0]
        p2 = nearest_edge[1]
        match_x, match_y = pedal((G.nodes[p1]['x'], G.nodes[p1]['y']), (G.nodes[p2]['x'], G.nodes[p2]['y']), (events.iloc[i]['x'], events.iloc[i]['y']))
        if(show_mod):
            print("Match ", match_x, match_y)
        
        #event_index = -1 * (max_index + i)
        edge_lengths = nx.get_edge_attributes(G, "length")
        originalLength = edge_lengths[(p1, p2, 0)]
        x_dif = G.nodes[p2]['x'] - G.nodes[p1]['x']
        y_dif = G.nodes[p2]['y'] - G.nodes[p1]['y']
        if not y_dif == 0:
            length_1e = (match_y - G.nodes[p1]['y'])/ (G.nodes[p2]['y'] - G.nodes[p1]['y']) * originalLength
        elif not x_dif== 0:
            length_1e = (match_x - G.nodes[p1]['x'])/ (G.nodes[p2]['x'] - G.nodes[p1]['x']) * originalLength
        else:
            length_1e = 0
        if length_1e < 1e-7:
            length_1e = 0
        length_2e = originalLength - length_1e
        if length_2e < 1e-7:
            length_2e = 0
        print(str(event_index) , " {" ,str(p1) , "： " ,str(length_1e) , ", " , str(p2) + ": " + str(length_2e) + "}")
        matched_points.append(event_index)
    return matched_points
        
def ntpi(t, pi, G,matched_index):
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=t, weight = 'length') #No need to start from scratch every time
    for p in paths.keys():
        if p in matched_index and (p not in found):
            #print(p)
            #print(p in found)
            found.append(p)
    
    return len(found)
def ntpi_grouped(t, pi, G,matched_index):
    count = 0
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=t, weight = 'length') #No need to start from scratch every time
    for p in paths.keys():
        if p in matched_index and (p not in found):
            #print(p)
            #print(p in found)
            count += int(G.nodes[p]['event_count'])
    
    return count
def Ltp (root, tj):
    D = 0
    #root = extended_shortest_path_tree(pi, G)
    for child in root.children:
        if  child.d <= tj:
            D = D + child.l + Ltp(child, tj)
        else:
            D = D + child.l + tj - child.d #if a node is further than tj, the subtree rooted at the node will not be computed
    return D
""" def wholeL():
    place = "Los Angeles, California"
    G = ox.graph_from_place(place, network_type='drive')
    fix_oneway(G)
    edge_lengths = nx.get_edge_attributes(G, "length")
    edge_oneway = nx.get_edge_attributes(G, "oneway")
    L = 0
    for edge in G.edges:
        L = L + edge_lengths[edge]/2
    return L """
def L(G):
    edge_lengths = nx.get_edge_attributes(G, "length")
    edge_oneway = nx.get_edge_attributes(G, "oneway")
    L = 0
    for edge in G.edges:
        L = L + edge_lengths[edge]/2
    return L
 
def Ktp(G, t, pi, n, matched_index):
    if n == 0:
        n = 1
    return L(G)/(n)*ntpi(t, pi, G, matched_index)
def Ktp_grouped(G, t, pi, n, matched_index):
    if n == 0:
        n = 1
    return L(G)/(n)*ntpi_grouped(t, pi, G, matched_index)
def Var(G, t, pi, n, et):
    l = L(G)
    Ltpi = Ltp(et, t)
    if n == 0:
        n = 1
    var = 1/(n) * l * Ltpi * (1-Ltpi/l)
    return var
def extended_shortest_path_tree(root_node, G, cutoff = None):
    root = Node(root_node, d = 0, l = 0) 
    edge_lengths = nx.get_edge_attributes(G, "length")
    if cutoff != None:
        cutoff += max(edge_lengths.values())
    distances, paths = nx.single_source_dijkstra(G,root_node,weight = 'length', cutoff = cutoff) #No need to restart 
    for pkey in paths.keys():
        if not pkey == root:
            path = paths[pkey]
            for i in range(1, len(path)):
                if not find(root, lambda node: node.name == path[i]):
                    if ( path[i-1], path[i], 0) in edge_lengths.keys():
                        node = Node(path[i], parent = find(root, lambda node: node.name == path[i-1]), l = edge_lengths[( path[i-1], path[i], 0)], d = distances[path[i]])
                    else:
                        node = Node(path[i], parent = find(root, lambda node: node.name == path[i-1]), l = edge_lengths[( path[i-1], path[i], 1)], d = distances[path[i]])

    for l in PreOrderIter(root):#Need to further confirm the correcteness
        lname = l.name
        neighborEdges = list(G.edges(lname))
        for e in neighborEdges:
            v2 = find(root, lambda node: node.name == e[1])
            if v2 and not (v2 == l.parent or l == v2.parent) :
                if (lname, e[1], 0) in edge_lengths.keys():
                    ledge = edge_lengths[(lname, e[1], 0)]
                else:
                    ledge = edge_lengths[(lname, e[1], 1)]
                #print('Leave: ', lname)
                #print('d: ', l.d)
                #print('neighbor d: ', v2.d)
                #print('Edge length: ', ledge)
                dbreak = (v2.d + ledge - l.d)/2
                node = Node(str(lname)+"break"+str(v2.name), parent = l, l = dbreak, d = l.d + dbreak)
                node = Node(str(v2.name)+"break"+str(lname), parent =v2, l = ledge- dbreak, d = v2.d + ledge- dbreak)
            #if not v2:
                #print("Error: Neighbor node not in tree: ",e, edge_lengths[(lname, e[1], 0)])
    return root

def extended_shortest_path_tree(root_node, G, cutoff = None):
    root = Node(root_node, d = 0, l = 0) 
    edge_lengths = nx.get_edge_attributes(G, "length")
    if cutoff != None:
        cutoff += max(edge_lengths.values())
    distances, paths = nx.single_source_dijkstra(G,root_node,weight = 'length', cutoff = cutoff) #No need to restart 
    for pkey in paths.keys():
        if not pkey == root:
            path = paths[pkey]
            for i in range(1, len(path)):
                if not find(root, lambda node: node.name == path[i]):
                    if ( path[i-1], path[i], 0) in edge_lengths.keys():
                        node = Node(path[i], parent = find(root, lambda node: node.name == path[i-1]), l = edge_lengths[( path[i-1], path[i], 0)], d = distances[path[i]])
                    else:
                        node = Node(path[i], parent = find(root, lambda node: node.name == path[i-1]), l = edge_lengths[( path[i-1], path[i], 1)], d = distances[path[i]])

    for l in PreOrderIter(root):#Need to further confirm the correcteness
        lname = l.name
        neighborEdges = list(G.edges(lname))
        for e in neighborEdges:
            v2 = find(root, lambda node: node.name == e[1])
            if v2 and not (v2 == l.parent or l == v2.parent) :
                if (lname, e[1], 0) in edge_lengths.keys():
                    ledge = edge_lengths[(lname, e[1], 0)]
                else:
                    ledge = edge_lengths[(lname, e[1], 1)]
                #print('Leave: ', lname)
                #print('d: ', l.d)
                #print('neighbor d: ', v2.d)
                #print('Edge length: ', ledge)
                dbreak = (v2.d + ledge - l.d)/2
                node = Node(str(lname)+"break"+str(v2.name), parent = l, l = dbreak, d = l.d + dbreak)
                node = Node(str(v2.name)+"break"+str(lname), parent =v2, l = ledge- dbreak, d = v2.d + ledge- dbreak)
            #if not v2:
                #print("Error: Neighbor node not in tree: ",e, edge_lengths[(lname, e[1], 0)])
    return root

def extended_shortest_path_tree_dict(root_node, G, cutoff = None):
    root = Node(root_node, d = 0, l = 0) 
    nodeDic = {}
    nodeDic[root_node] = root;
    edge_lengths = nx.get_edge_attributes(G, "length")
    if cutoff != None:
        cutoff += max(edge_lengths.values())
    distances, paths = nx.single_source_dijkstra(G,root_node,weight = 'length', cutoff = cutoff) #No need to restart 
    
    for pkey in paths.keys():
        if not pkey == root:
            path = paths[pkey]
            for i in range(1, len(path)):
                if not path[i] in nodeDic:
                    if ( path[i-1], path[i], 0) in edge_lengths.keys():
                        node = Node(path[i], parent = nodeDic[path[i-1]], l = edge_lengths[( path[i-1], path[i], 0)], d = distances[path[i]])
                    else:
                        node = Node(path[i], parent = nodeDic[path[i-1]], l = edge_lengths[( path[i-1], path[i], 1)], d = distances[path[i]])
                    nodeDic[path[i]] = node
    for l in PreOrderIter(root):#Need to further confirm the correcteness
        lname = l.name
        neighborEdges = list(G.edges(lname))
        for e in neighborEdges:
            if e[1] in nodeDic.keys():
                v2 = nodeDic[e[1]]
                if not (v2 == l.parent or l == v2.parent) :
                    if (lname, e[1], 0) in edge_lengths.keys():
                        ledge = edge_lengths[(lname, e[1], 0)]
                    else:
                        ledge = edge_lengths[(lname, e[1], 1)]
                    dbreak = (v2.d + ledge - l.d)/2
                    node = Node(str(lname)+"break"+str(v2.name), parent = l, l = dbreak, d = l.d + dbreak)
                    node = Node(str(v2.name)+"break"+str(lname), parent =v2, l = ledge- dbreak, d = v2.d + ledge- dbreak)
    return root

def remove_node(G, node):
    out_edges = list(G.edges(node))
    if(len(out_edges) != 2):
        print("Error: Number of ut edges for node ", str(node), " is not 2")
        return
    node1 = out_edges[0][1]
    node2 = out_edges[1][1]
    length1 = G.edges[(node, node1, 0)]['length']
    length2 = G.edges[(node, node2, 0)]['length']
    rec_length = length1 + length2
    G.remove_node(node)
    G.add_edges_from(((node1, node2, 0), (node2,node1, 0)))
    attrs = {(node1, node2, 0):{"length": rec_length, "oneway": False}, (node2, node1, 0):{"length": rec_length, "oneway": False}}
    nx.set_edge_attributes(G, attrs)
def shortest_path_tree(root_node, G):
    root = Node(root_node, d = 0, l = 0)
    edge_lengths = nx.get_edge_attributes(G, "length")
    distances, paths = nx.single_source_dijkstra(G,root_node,weight = 'length') #No need to start 
    for pkey in paths.keys():
        if not pkey == root:
            path = paths[pkey]
            for i in range(1, len(path)):
                if not find(root, lambda node: node.name == path[i]):
                    node = Node(path[i], parent = find(root, lambda node: node.name == path[i-1]), l = edge_lengths[( path[i-1], path[i], 0)], d = distances[path[i]])
    return root

def Var(G, t, pi, n, et):
    l = L(G)
    Ltpi = Ltp(et, t)
    var = 1/(n-1) * l * Ltpi * (1-Ltpi/l)
    return var
def Var_(G, t, pi, n, et):
    l = L(G)
    Ltpi = Ltp(et, t)
    var = 1/(n-1) * l * Ltpi * (1-Ltpi/l)
    return var
def detectHotspot(index, crimes, t = 10, m = 20, plot = True, show_log = False):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    p = crimes.iloc[index]
    location_point = (p['LAT'], p['LON'])
    wn = (p['LON'] - 0.03, p['LAT'] + 0.03 )
    es = (p['LON'] + 0.03, p['LAT'] - 0.03 )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=ox.graph_from_bbox(wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        return 0
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['LAT'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LAT'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['LON'], 'y': crimes_in_range['LAT']})

    center_index = -1
    for i in range(0, crimes_in_range.shape[0]):
        if crimes_in_range.iloc[i]['LAT'] == p['LAT'] and crimes_in_range.iloc[i]['LON'] == p['LON']:
            center_index = i
            break

    matched_index = match_event_to_edges(events, G)
    pi = matched_index[center_index]
    if plot:
        nc = ["r" if (node in matched_index and node != matched_index[center_index]) else "b" for node in G.nodes()]
        ns = [10 if (node in matched_index) else 0 for node in G.nodes()]
        fig, ax = ox.plot_graph(G, node_color=nc, node_size = ns, bgcolor = '#ffffff', save = True, filepath = "./" + str(c_code) + "_" + str(t) + "_" + str(index) + ".jpg")
    time1 = time.time()
    et = extended_shortest_path_tree(pi, G, cutoff = t * m) 
    #pi = matched_index[center_index]
    cluster = 0
    file1 = open("cluster_northLA.txt", "a")  # append mode
    file1.write(str(index)+ ", ")
    for i in range(1, m):
        print("i = ", i)
        stime = time.time()
        mean = Ltp(et, t * i)
        ltime = time.time()
        print("Ltp cal: ", mean)
        std = Var(G, t * i, pi, len(matched_index)-1, et)**0.5
        sTtime = time.time()
        print("std cal: ", std)
        k = Ktp(G, t * i, pi, len(matched_index)-1, matched_index)
        ktime = time.time()
        print("k cal: ",k)
        normalDist = scipy.stats.norm(mean, std)
        if k > normalDist.ppf(1-.01):
            file1.write("1, ")
            print("At ti = ", t*i, " k >upper critical")
            cluster = cluster + 1
        if k < normalDist.ppf(.01):
            file1.write("-1, ")
            print("At ti = ", t*i, " k <lower critical")
            cluster = cluster - 1
    time2 = time.time()
    
    file1.write(str(time2 - time1) + "\n")
    file1.close()
    return G
def detectHotspot_prob(prob, crimes, t = 10, m = 20, plot = True, show_log = False):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    #p = crimes.iloc[index]
    #location_point = (p['LAT'], p['LON'])
    print(prob)
    wn = (prob[0] - 0.03, prob[1] + 0.03 )
    es = (prob[0] + 0.03, prob[1] - 0.03 )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=ox.graph_from_bbox(wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        return 0
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['LAT'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LAT'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['LON'], 'y': crimes_in_range['LAT']})
    events = events.append({'x': prob[0], 'y': prob[1]}, ignore_index = True)
    if(show_log):
        print("No. of events: ", len(events))
        print(events)
    #center_index = events.shape[0] - 1
    #center_index = -1
    #for i in range(0, crimes_in_range.shape[0]):
        #if crimes_in_range.iloc[i]['LAT'] == p['LAT'] and crimes_in_range.iloc[i]['LON'] == p['LON']:
            #center_index = i
            #break

    matched_index = match_event_to_edges(events, G, show_mod = show_log)
    if(show_log):
        print("No. of matched nodes: ", len(matched_index))
    center_index = len(matched_index)-1
    pi = matched_index[center_index]
    time1 = time.time()
    et = extended_shortest_path_tree(pi, G, cutoff = t * m) 
    #pi = matched_index[center_index]
    cluster = 0
    file1 = open("cluster_midLA_-118.5600, 33.9253.txt", "a")  # append mode
    file1.write(str(prob)+ ", ")
    if(len(events) == 1):
        if(show_log):
            print("No event")
        file1.write("No event\n")
        return
    for i in range(1, m):
        print("i = ", i)
        print("No. of points in range: ", ntpi(t*i, pi, G,matched_index))
        stime = time.time()
        mean = Ltp(et, t * i)
        ltime = time.time()
        print("Ltp cal: ", mean)
        std = Var(G, t * i, pi, len(matched_index), et)**0.5
        sTtime = time.time()
        print("std cal: ", std)
        k = Ktp(G, t * i, pi, len(matched_index), matched_index)
        ktime = time.time()
        print("k cal: ",k)
        normalDist = scipy.stats.norm(mean, std)
        if k > normalDist.ppf(1-.01):
            file1.write("1, ")
            print("At ti = ", t*i, " k >upper critical")
            cluster = cluster + 1
        elif k < normalDist.ppf(.01):
            file1.write("-1, ")
            print("At ti = ", t*i, " k <lower critical")
            cluster = cluster - 1
        else:
            file1.write("0, ")
            print("At ti = ", t*i, " CSR")
            #cluster = cluster - 1
    time2 = time.time()
    if plot:
        nc = ["r" if (node in matched_index and node != matched_index[center_index]) else "b" for node in G.nodes()]
        ns = [10 if (node in matched_index) else 0 for node in G.nodes()]
        classPath = str(t) + "_" + str(m)+"/"
        if cluster > 0:
            classPath  = classPath + "Cluster"
        elif cluster < 0:
            classPath = classPath + "De_cluster"
        else:
            classPath = classPath + "CSR"
        fig, ax = ox.plot_graph(G, node_color=nc, node_size = ns, bgcolor = '#ffffff', save = True, filepath = "./" + classPath + "/" + str(c_code) + "_" + str(t) + "_" + str(prob) + ".jpg")
    file1.write(str(cluster) + ", "+ str(time2 - time1) + "\n")
    file1.close()
    remove_node(G, center_index)
    return G

def map_events_to_tile_grouped(prob, crimes,  base_id = 0, show_log = False):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    #p = crimes.iloc[index]
    #location_point = (p['LAT'], p['LON'])
    #print(prob)
    wn = (prob[0] - 0.031, prob[1] + 0.031 )
    es = (prob[0] + 0.031, prob[1] - 0.031 )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=ox.graph_from_bbox(wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        raise
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['Latitude'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['Latitude'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['Longitude'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['Longitude'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['Longitude'], 'y': crimes_in_range['Latitude']})
    #events = events.append({'x': prob[0], 'y': prob[1]}, ignore_index = True)
    if(show_log):
        print("No. of events: ", len(events))
        print(events)
    matched_index = match_event_to_edges_grouped(events, G, show_mod = show_log)
    return G
def map_events_to_tile(prob, crimes,  base_id = 0, show_log = False):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    #p = crimes.iloc[index]
    #location_point = (p['LAT'], p['LON'])
    #print(prob)
    wn = (prob[0] - 0.031, prob[1] + 0.031 )
    es = (prob[0] + 0.031, prob[1] - 0.031 )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=ox.graph_from_bbox(wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        return 0
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['LAT'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LAT'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['LON'], 'y': crimes_in_range['LAT']})
    #events = events.append({'x': prob[0], 'y': prob[1]}, ignore_index = True)
    if(show_log):
        print("No. of events: ", len(events))
        print(events)
    matched_index = match_event_to_edges(events, G, show_mod = show_log)
    return G

def map_events_to_tile_distance(prob, crimes,  base_id = 0, show_log = False):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    #p = crimes.iloc[index]
    #location_point = (p['LAT'], p['LON'])
    #print(prob)
    wn = (prob[0] - 0.031, prob[1] + 0.031 )
    es = (prob[0] + 0.031, prob[1] - 0.031 )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=ox.graph_from_bbox(wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        return 0
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['LAT'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LAT'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['LON'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['LON'], 'y': crimes_in_range['LAT']})
    #events = events.append({'x': prob[0], 'y': prob[1]}, ignore_index = True)
    if(show_log):
        print("No. of events: ", len(events))
        print(events)
    matched_index = match_event_to_edges_distance(events, G, show_mod = show_log)
    return G

def network_distance_count_node(pi, max_dist, step_size, G, grouped):
    n_list = [0] * (max_dist//step_size)
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=max_dist, weight = 'length') #No need to start from scratch every time
    for p in distances.keys():
        if p < 0 and (p not in found):
            #print(p)
            #print(p in found)
            if grouped:
                n_list[int(distances[p]//step_size)] += G.nodes[p]['event_count']
            else:
                n_list[int(distances[p]//step_size)] += 1
            
    return n_list


def network_distance_count_pos(center_point, max_dist, step_size, G, grouped = False):
    p = ox.distance.nearest_nodes(G, center_point[0], center_point[1])
    return network_distance_count_node(p, max_dist, step_size, G, grouped)

def fix_negative_edge(G):
    edge_lengths = nx.get_edge_attributes(G, "length")
    edge_lengths
    for edge in edge_lengths.keys():
        if edge_lengths[edge] < 0:
            attrs = attrs = {edge:{"length": 0, "oneway": False}}
            nx.set_edge_attributes(G, attrs)

def wholeL(G):
    edge_lengths = nx.get_edge_attributes(G, "length")
    L = 0
    for edge in G.edges:
        L = L + edge_lengths[edge]/2
    return L
def wholeN(G):
    count = 0
    for node in G.nodes:
        if node < 0:
            count += 1
    return count
    
def detectHotspot_matched_prob(prob, G, crimes, t = 20, m = 80, ccode = "default", max_distance = 4800, global_HS = False, plot = True, show_log = False):
    c_code = "All"
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
def detectHotspot_matched_prob(prob, G, t = 20, m = 80, c_code = "default", max_distance = 4800, global_HS = False, plot = True, show_log = False):
    #c_code = "All"
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
def euclidean_distance(G, u, v):
    x1 = G.nodes[u]["x"]
    x2 = G.nodes[v]["x"]
    y1 = G.nodes[u]["y"]
    y2 = G.nodes[v]["y"]
    return geopy.distance.geodesic((y1, x1), (y2, x2)).m

def fix_edge_error_after_merging_grids(G):
    length = nx.get_edge_attributes(G, "length")
    error_edge_list = list()
    error_count = 0
    for u, v, k in G.edges(keys=True):
        if (u < 0 or v < 0) and abs(euclidean_distance(G, u,v) - length[u, v, k]) > 1:
            error_count += 1
            error_edge_list.append((u,v,k))
    G.remove_edges_from(error_edge_list)
    print(str(error_count), " edges removed")

def compute_accuracy(y_true, y_pred):
    correct_predictions = 0
    # iterate over each label and check
    for true, predicted in zip(y_true, y_pred):
        if true == predicted:
            correct_predictions += 1
    # compute the accuracy
    accuracy = correct_predictions/len(y_true)
    return accuracy
def network_distance_count_lengthProp(pi, max_dist, step_size, G):
    n_list = [0] * (max_dist//step_size)
    n_list_prop = [0] * (max_dist//step_size)
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=max_dist, weight = 'length') #No need to start from scratch every time
    for p in distances.keys():
        if p < 0 and (p not in found):
            #print(p)
            #print(p in found)
            n_list[int(distances[p]//step_size)] += 1
    for i in range(0, (max_dist//step_size)):
        n_list_prop[i] = n_list[i] / ((2 * (i+1) -1))
    return n_list_prop

def network_distance_count_lengthPropSimp(pi, max_dist, step_size, G):
    n_list = [0] * (max_dist//step_size)
    n_list_prop = [0] * (max_dist//step_size)
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=max_dist, weight = 'length') #No need to start from scratch every time
    for p in distances.keys():
        if p < 0 and (p not in found):
            #print(p)
            #print(p in found)
            n_list[int(distances[p]//step_size)] += 1
    for i in range(0, (max_dist//step_size)):
        n_list_prop[i] = n_list[i] / ((i+1))
    return n_list_prop

def network_distance_count_posProp(center_point, max_dist, step_size, G):
    p = ox.distance.nearest_nodes(G, center_point[0], center_point[1])
    return network_distance_count_lengthProp(p, max_dist, step_size, G)

def network_distance_count_posPropSimp(center_point, max_dist, step_size, G):
    p = ox.distance.nearest_nodes(G, center_point[0], center_point[1])
    return network_distance_count_lengthPropSimp(p, max_dist, step_size, G)

def network_distance_countRoadProp(center_point, max_dist, step_size, G):
    pi = ox.distance.nearest_nodes(G, center_point[0], center_point[1])
    n_list = [0] * (max_dist//step_size)
    n_list_prop = [0] * (max_dist//step_size)
    found = []
    distances, paths = nx.single_source_dijkstra(G,pi,cutoff=max_dist, weight = 'length') #No need to start from scratch every time
    for p in distances.keys():
        if p < 0 and (p not in found):
            #print(p)
            #print(p in found)
            n_list[int(distances[p]//step_size)] += 1

    dist_to_prev = {}
    dist_list = [0] * (max_dist//step_size)
    keys = list(paths)
    """for key in keys:
        if(key > 0 and key != pi):
            del paths[key]
            del distances[key]
        else:
            for node in paths[key]:
                if (node > 0 and node != pi):
                    paths[key].remove(node)"""
    for key,value in paths.items():
        #value = [x for x in value if x < 0 or x == pi]
        if(len(value) <= 1):
            dist_to_prev[key] = distances[key]
        else:
            #print(value)
            parent_node = value[-2]
            dist_to_prev[key] = distances[key] - distances[parent_node]
            pos_child = int(distances[key] / step_size)
            pos_parent = int(distances[parent_node] / step_size)
            ceil_parent = math.ceil(distances[parent_node] / step_size)
            if(ceil_parent < pos_child):
                for i in range(ceil_parent, pos_child-1):
                    dist_list[i] = dist_list[i] + step_size
            if pos_child > pos_parent:
                dist_list[pos_parent] += (ceil_parent) * step_size - distances[parent_node]
            dist_list[pos_child] += distances[key] - pos_child * step_size
    for i in range(0, (max_dist//step_size)):
        if(dist_list[i] == 0 and n_list[i] == 0 ):
            n_list_prop[i] = 0
        elif dist_list[i] == 0 and n_list[i] > 0:
            print("Error: division by 0 at", pi)
            print(n_list[i], i)
            print(distances)
        else:
            n_list_prop[i] = n_list[i] / dist_list[i]
    return n_list_prop

def partition(G, no_of_steps = 4):
    nodes = ox.utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]
    events = nodes[nodes.index<0]
    #no_of_step = 4
    grids = []
    for i in range(0, no_of_step):
        row = []
        for i in range(0, no_of_step):
            row.append([])
        grids.append(row)
    for index, row in events.iterrows():
        #print(row['x'], row['y'])
        x = int((row['x'] - min_x)/step_x)
        y = int((row['y'] - min_y)/step_y)
        if(x == no_of_step):
            x = x-1
        if(y == no_of_step):
            y = y -1
        grids[x][y].append(index)
    return grids

def detectHotspot_matched_prob_gridValidation(prob, G, t = 20, m = 80, c_code = "default", max_distance = 4800, global_HS = False, plot = True, show_log = False):
    #c_code = "All"
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
    return cluster

def encode_data_list(G, data_dir, max_dist, step_size, show_log = False):
    encode_list= list()
    file_count = 0
    total_time = 0
    for filename in os.listdir(data_dir):
        file_count += 1
        words = filename.split("_")
        coord = words[3]
        res = tuple(map(float, coord[1:-5].split(",")))
        es_time = time.time()
        encode = network_distance_count_pos(res, max_dist, step_size, G)
        ee_time = time.time()
        total_time += ee_time - es_time
        if(show_log):
            print("Time for feature extraction for " + filename + " is " + str(ee_time - es_time))
        encode_list.append(encode)
    return encode_list, (total_time / file_count)

def encode_data_list_prop(G, data_dir, max_dist, step_size, show_log = False):
    encode_list= list()
    file_count = 0
    total_time = 0
    for filename in os.listdir(data_dir):
        file_count += 1
        words = filename.split("_")
        coord = words[3]
        res = tuple(map(float, coord[1:-5].split(",")))
        es_time = time.time()
        encode = network_distance_count_posProp(res, max_dist, step_size, G)
        ee_time = time.time()
        total_time += ee_time - es_time
        if(show_log):
            print("Time for feature extraction for " + filename + " is " + str(ee_time - es_time))
        encode_list.append(encode)
    return encode_list, (total_time / file_count)

def encode_data_list_propSimp(G, data_dir, max_dist, step_size, show_log = False):
    encode_list= list()
    file_count = 0
    total_time = 0
    for filename in os.listdir(data_dir):
        file_count += 1
        words = filename.split("_")
        coord = words[3]
        res = tuple(map(float, coord[1:-5].split(",")))
        es_time = time.time()
        encode = network_distance_count_posPropSimp(res, max_dist, step_size, G)
        ee_time = time.time()
        total_time += ee_time - es_time
        if(show_log):
            print("Time for feature extraction for " + filename + " is " + str(ee_time - es_time))
        encode_list.append(encode)
    return encode_list, (total_time / file_count)

def encode_data_list_grouped(G, data_dir, max_dist, step_size, show_log = False):
    encode_list= list()
    file_count = 0
    total_time = 0
    for filename in os.listdir(data_dir):
        file_count += 1
        words = filename.split("_")
        coord = words[3]
        res = tuple(map(float, coord[1:-5].split(",")))
        es_time = time.time()
        encode = network_distance_count_pos_grouped(res, max_dist, step_size, G)
        ee_time = time.time()
        total_time += ee_time - es_time
        if(show_log):
            print("Time for feature extraction for " + filename + " is " + str(ee_time - es_time))
        encode_list.append(encode)
    return encode_list, (total_time / file_count)
def map_events_to_tile_cropped_by_place(place, prob, crimes,  group = False, base_id = 0, show_log = False):
    wn = (prob[0] - 0.031, prob[1] + 0.031 )
    es = (prob[0] + 0.031, prob[1] - 0.031 )
    try:
        G=graph_from_bbox_with_place(place, wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=False)
    except Exception as e:
        raise
    #G=ox.graph_from_point(location_point, dist = 1609, network_type="drive_service", simplify=False)
    crimes_in_range = crimes[crimes['Latitude'] > es[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['Latitude'] < wn[1]]
    crimes_in_range = crimes_in_range[crimes_in_range['Longitude'] > wn[0] ]
    crimes_in_range = crimes_in_range[crimes_in_range['Longitude'] < es[0] ]
    #if show_log:
       # print("No. of crimes: ", crimes_in_range.shape[0])
        #print(crimes_121_in_range)
    fix_oneway(G)
    events = pd.DataFrame({'x': crimes_in_range['Longitude'], 'y': crimes_in_range['Latitude']})
    #events = events.append({'x': prob[0], 'y': prob[1]}, ignore_index = True)
    if(show_log):
        print("No. of events: ", len(events))
        print(events)
    if (group):
        matched_index = match_event_to_edges_grouped(events, G, show_mod = show_log)
    else:
        matched_index = match_event_to_edges(events, G, show_mod = show_log)
    return G
def graph_from_bbox_with_place(
    query,
    north=None,
    south=None,
    east=None,
    west=None,
    #bbox=None,
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=None,
    custom_filter=None,
):
    gdf_place = ox.geocoder.geocode_to_gdf(
                query, which_result=None, buffer_dist=None
            )
    place_polygon = gdf_place["geometry"].unary_union
    
    bbox = (north, south, east, west)
    bbox_polygon = ox.utils_geo.bbox_to_poly(bbox=bbox)

    intersection_polygon = shapely.intersection(place_polygon, bbox_polygon)
    G = ox.graph_from_polygon(
            intersection_polygon,
            network_type=network_type,
            simplify=simplify,
            retain_all=retain_all,
            truncate_by_edge=truncate_by_edge,
            clean_periphery=clean_periphery,
            custom_filter=custom_filter,
        )
    return G