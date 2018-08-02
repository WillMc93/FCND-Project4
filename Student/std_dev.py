import csv
import numpy as np

path = '/home/will/Desktop/FCND-Estimation-CPP/config/log/Graph1.txt'

data = list()

with open(path, 'r') as file:
	file.readline() # throw away labels

	datareader = csv.reader(file)
	for row in datareader:
		data.append(float(row[1]))
	data = np.array(data) # I forgot to cast to float in an earlier version \
						  # and believe this line to be superflous

print("Std_dev = {}".format(np.std(data)))
	