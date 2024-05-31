from shapely.geometry import Point
import multiprocessing
from hotspot import *
from os.path import exists
import sys
import os
show_log = False
def report_status(result):
    if(show_log):
        print(f'Callback received: {result}')
def map_and_save(place, folder, prob, g_id, crimes):
    #folder = "GraphML_CH_30days_cropped"
    #place = "Chicago"
    if(exists(folder + "/" + str(prob) + ".graphml")):
        if(show_log):
            print("Already exists, skip " + str(prob))
        return 'skip'
    else:
        if(show_log):
            print("Create area centered at ", str(prob))
    try:
        G = match_events_to_tile_cropped_by_place(place, prob, crimes, base_id = crimes.shape[0] * g_id, show_log = True)
        #print(result)
        #GList.append(G)
        ox.io.save_graphml(G, filepath = folder + "/" + str(prob) + ".graphml")
    except Exception as e:
        return "Fail to match points in grid " + str(prob) + ": " + str(e)
    return "Grid centered at " + str(prob) + " created successfully!"
def check_existence(folder, prob, g_id, crimes):
    if(exists(folder + "/" + str(prob) + ".graphml")):
        print("Already exists, skip " + str(prob))
    else:
        print(str(prob) + "does not exists")
def match_points_to_network(location, data_file, lon, lat, number_of_processes, output_folder = "Temp_Grids"):
    points = pd.read_csv(data_file)
    if(show_log):
        print("Total number of events: " + str(points.shape[0]))
    xmin = points[lon].min()
    xmax = points[lon].max()
    ymin = points[lat].min()
    ymax = points[lat].max()
    p = multiprocessing.Pool(int(number_of_processes))
    #prob = (-118.3621635,34.0122962)
    #GList = list()
    i = 100
    x = xmin
    while x < xmax:
        y = ymin
        while y < ymax:
            p.apply_async(map_and_save, (location, output_folder, (x,y), i, points,), callback=report_status)
            #map_and_save((x,y), i, crimes)
            #message = poolResult.get()
            #print(message)
            #print(map_and_save((x,y), i, crimes))
            y = y + 0.06
            i = i+1
        x = x + 0.06
    p.close()
    p.join() 
    return combine_grid(output_folder)
def combine_grid(graphml_dir):
    GridList = list()
    for filename in os.listdir(graphml_dir):
        Gload = ox.io.load_graphml(filepath = graphml_dir + "/"+filename)
        GridList.append(Gload)
    GridM =  nx.compose_all(GridList)
    fix_edge_error_after_merging_grids(GridM)
    fix_negative_edge(GridM)
    return GridM
if __name__ == '__main__': 
    location = sys.argv[1]
    data_file = sys.argv[2]
    output_folder = sys.argv[3]
    lon = sys.argv[4]
    lat = sys.argv[5]
    np = sys.argv[6]
    
def match_events_to_tile_cropped_by_place(place, prob, crimes, insert = False, group = False, base_id = 0, show_log = False, tail_size = 0.031):
    #wn = (-118.2730593,34.0532872)
    #es = (-118.2121635,34.0122962)
    #p = crimes.iloc[index]
    #location_point = (p['LAT'], p['LON'])
    #print(prob)
    wn = (prob[0] - tail_size, prob[1] + tail_size )
    es = (prob[0] + tail_size, prob[1] - tail_size )
    #location_point = (34.0463279, -118.2678461)
    try:
        G=graph_from_bbox_with_place(place, wn[1], es[1], es[0], wn[0], network_type="drive_service", simplify=True)
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
        matched_index = match_events_to_edges_grouped(events, G, insert = insert, show_mod = show_log)
    else:
        matched_index = match_events_to_edges(events, G, insert = insert, show_mod = show_log)
    return G

def match_events_to_edges(events, G, group = False,insert = False, show_mod = False): #May need to create a function for reverting the ops
    matched_points = []
    for i in range(len(events)):
        #event_index = -1 * (max_index + i)
        event_index = events.iloc[i].name * -1
        if (show_mod):
            print("Finding nearest edge for ", event_index, ", ", events.iloc[i])
        nearest_edge = ox.nearest_edges(G, events.iloc[i]['x'], events.iloc[i]['y'])
        if (show_mod):
            print("Nearest edge is ", nearest_edge)
        #Compute the projected position
        if('geometry' in G.edges[nearest_edge]):#Simplified graph
            point_geom = Point(events.iloc[i]['y'], events.iloc[i]['x'], )
            edge_geom = G.edges[nearest_edge]['geometry']
            nearest_point_on_edge = edge_geom.interpolate(edge_geom.project(point_geom))
            #print(nearest_point_on_edge)
            match_x, match_y = nearest_point_on_edge.x, nearest_point_on_edge.y
            if(show_mod):
                    print("Matched through interpolate ", match_x, match_y)
        else:#Straigntline or unsimplified graph
            p1 = nearest_edge[0]
            p2 = nearest_edge[1]
            match_x, match_y = pedal((G.nodes[p1]['x'], G.nodes[p1]['y']), (G.nodes[p2]['x'], G.nodes[p2]['y']), (events.iloc[i]['x'], events.iloc[i]['y']))
            if(show_mod):
                print("Matched through pedal ", match_x, match_y)

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
        if(insert):
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
        else:
            if('event_count' in G.edges[nearest_edge]):
                G.edges[nearest_edge]['event_count'] += 1
            else:
                G.edges[nearest_edge]['event_count'] = 1
            distance_low = -1
            if p1 < p2:
                distance_low = length_1e
            else:
                distance_low = length_2e
            if 'event_dist_to_smaller_node' in G.edges[nearest_edge]:
                G.edges[nearest_edge]['event_dist_to_smaller_node'][event_index] = distance_low
            else:
                G.edges[nearest_edge]['event_dist_to_smaller_node'] = {event_index: distance_low}
            if 'event_projected_coordinate' in G.edges[nearest_edge]:
                G.edges[nearest_edge]['event_projected_coordinate'][event_index] = (match_x, match_y)
            else:
                G.edges[nearest_edge]['event_projected_coordinate']= {event_index: (match_x, match_y)}
            
        matched_points.append(event_index)
    return matched_points

def match_event_to_edges_grouped(events, G, group_dist = 10, show_mod = False): #May need to create a function for reverting the ops
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