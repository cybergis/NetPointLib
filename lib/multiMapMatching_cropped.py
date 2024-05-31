#Documents\Hotspot\network-clustering-main\dataprepare\Templates
#python multiMapMatching_cropped.py [location] [datafile]

import multiprocessing
from hotspot import *
from os.path import exists
import sys

def report_status(result):
	print(f'Callback received: {result}')
def map_and_save(place, folder, prob, g_id, crimes):
	#folder = "GraphML_CH_30days_cropped"
	#place = "Chicago"
	if(exists(folder + "/" + str(prob) + ".graphml")):
		print("Already exists, skip " + str(prob))
		return 'skip'
	else:
		print("Create area centered at ", str(prob))
	try:
		G = map_events_to_tile_cropped_by_place(place, prob, crimes, base_id = crimes.shape[0] * g_id, show_log = False)
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
if __name__ == '__main__': 
	location = sys.argv[1]
	data_file = sys.argv[2]
	output_folder = sys.argv[3]
	lon = sys.argv[4]
	lat = sys.argv[5]
	np = sys.argv[6]
	points = pd.read_csv(data_file)
	print("Total number of events: " + str(points.shape[0]))
	xs = points[lon]
	ys = points[lat]
	xmin = points[lon].min()
	xmax = points[lon].max()
	ymin = points[lat].min()
	ymax = points[lat].max()
	p = multiprocessing.Pool(int(np))
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
	print("Done!")
	#GM = nx.compose_all(GList)
	#ox.io.save_graphml(GM, filepath = "GraphML/Merged.graphml")
