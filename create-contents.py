import os
import numpy as np
import json
import sys

#search following directory
indir = '/home/filip/Documents/www.ngdc.noaa.gov/ns41'

#create new file containing summary of the contents
contents = open("contenst_ns41.txt" , "w+")
contents.write("file_name   year   min_decimal_day   max_decimal_day   min_Geographic_Latitude   max_Geographic_Latitude   min_Geographic_Longitude   max_Geographic_Longitude   min_Rad_Re   max_Rad_Re   \n")

#opens each file in the directory skipping first 5 lines
for root, dirs, filenames in os.walk(indir):
   for plik in filenames:
      log = open(os.path.join(root, plik), 'r')
      reading = log.readlines()[5:] 

#creats a file header to store info about columns for one iteration
      header = open("head.json", "w+")
      for line in reading:
         if line.startswith('#'):
            li = line.strip('#')
            header.write(li)
      header.close()
      
      head = json.load(open("head.json")) #load header file

      f = np.loadtxt(os.path.join(root, plik)) #use numpy library to create arrays

      variables = {}
      for x, value in head.items():
         i = value['START_COLUMN']
         j = i + value['DIMENSION'][0]
         variables.update({x: f[ : , i : j] for column in head})

#example of the data we might want to extract from the folder to compare with other satellites
      min_decimal_day = np.min(variables["decimal_day"])
      max_decimal_day = np.max(variables["decimal_day"])
      min_Geographic_Latitude = np.min(variables["Geographic_Latitude"])
      max_Geographic_Latitude = np.max(variables["Geographic_Latitude"])
      min_Geographic_Longitude = np.min(variables["Geographic_Longitude"])
      max_Geographic_Longitude = np.max(variables["Geographic_Longitude"])
      min_Rad_Re = np.min(variables["Rad_Re"])
      max_Rad_Re = np.max(variables["Rad_Re"])
      year = np.max(variables["year"])
      
#save the info to the file
      contents.write("%s   %f   %f   %f   %f   %f   %f   %f   %f   %f   \n" % (plik , year , min_decimal_day , max_decimal_day , min_Geographic_Latitude , max_Geographic_Latitude , min_Geographic_Longitude , max_Geographic_Longitude , min_Rad_Re , max_Rad_Re ))

contents.close()


