import numpy as np
import json
from collections import OrderedDict

fil = open("ns57_080113_v1.03.ascii")
reading = fil.readlines()[5:] 

header = open("head.json", "w+")
for line in reading:
   if line.startswith('#'):
      li = line.strip('#')
      header.write(li)

header.close()

#Use the header file to get START_COLUMN, DIMENSION and then slice the data into arrays.

head = json.load(open("head.json") , object_pairs_hook = OrderedDict)

f = np.loadtxt("ns57_080113_v1.03.ascii")

variables = {}
for x, value in head.items():
   i = value['START_COLUMN']
   j = i + value['DIMENSION'][0]
   variables.update({x: f[ : , i : j] for column in head})

print variables["pfitpars"]
