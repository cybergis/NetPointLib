#Documents\Hotspot\network-clustering-main\dataprepare\Templates
import multiprocessing
from hotspot import *
from os.path import exists
def report_status(result):
	print(f'Callback received: {result}')
def map_and_save(prob, g_id, crimes):
	folder = "GraphML_CH_30days_cropped"
	place = "Chicago"
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
def check_existence(prob, g_id, crimes):
	if(exists("GraphML50k_new/" + str(prob) + ".graphml")):
		print("Already exists, skip " + str(prob))
	else:
		print(str(prob) + "does not exists")
if __name__ == '__main__': 
	crimesAll = pd.read_csv("Crimes_-_last 30.csv")
	#c_code = 510
	#crimes = crimesAll[(crimesAll['Crm Cd'] == 330) | (crimesAll['Crm Cd'] == 740)] 
	crimes = crimesAll
	print("Total number of events: " + str(crimes.shape[0]))
	xs = crimes['Longitude']
	ys = crimes['Latitude']
	p = multiprocessing.Pool(4)
	#prob = (-118.3621635,34.0122962)
	#GList = list()
	i = 100
	x = -87.9074726 
	while x < -87.52461633:
		y = 41.64458971
		while y < 42.02253615:
			p.apply_async(map_and_save, ((x,y), i, crimes,), callback=report_status)
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
