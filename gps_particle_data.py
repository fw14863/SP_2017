# For me (Ubuntu) I needed to install the wget module (details are here: https://ubuntuforums.org/showthread.php?t=2287729) 
import wget
import os
import numpy as np 
import time
import json
import sys
import math
import re
from io import StringIO

from datetime import datetime, timedelta
from itertools import compress

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from urllib2 import urlopen

from mpl_toolkits.basemap import Basemap

import aacgmv2


def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(days=n*0.25);

#====================================================================================================
#
# This class holds all of the data and functions associated with a single data file
#
#
#====================================================================================================
class gps_particle_datafile():
    # Variables in the data class 
    
    def __init__(self,filepath, sat_number):
        # Generate the basepath and filename 
        self.satellite_number = sat_number;
        # self.basepath = self.folder_path();
        self.data_filename = filepath;
        self.data = self.webread();        
        
        # Create the JSON header file
        # self.json_filename = self.create_json();
             

#==============================================================================
#     def __init__(self,sat_number,day,month,year):
#         
#         self.satellite_number = sat_number;
#         self.day = day;
#         self.month = month;
# 
#         if(year>1990):
#             self.date = datetime(year,month,day);
#             self.year = year-2000;
#         
#         if((year>20) and (year<1990)):
#             self.date = datetime(year+1900,month,day);
#             self.year = year;
#         elif((year<20) and (year<1990)):
#             self.date = datetime(year+2000,month,day);
#             self.year = year;
#             
#         # Generate the basepath and filename 
#         self.basepath = self.folder_path();
#         self.data_filename = self.file_path();
#         self.fail_load=False;
# 
#         # Check if the file is available locally in this folder?
#         exist_test = os.path.isfile(self.data_filename);
#         if(not exist_test):
#             # No it's been downloaded so download it
#             print ' - Downloading ' + self.data_filename;
#             if(self.download()):
#                 self.fail_load=True;
#                 return;
#         else:
#             print ' - File ' + self.data_filename + ' already exists ';
#             
#         if(not self.fail_load):
#     
#             # Create the JSON header file
#             self.json_filename = self.create_json();
#             
#             # Read the data from the file 
#             self.data = self.read();
#             
#             # self.summary.print_summary();
#             # print self.get_next_datecode();
#             # print self.get_last_date();
#             
#==============================================================================
    # Returns the end time of the data file
    def get_last_date(self):

        # Get the decimal days
        doys = self.data['decimal_day'];
        
        # Check if the file goes over a new year boundary 
        new_year = any(doys>365);
        
        if(new_year): 
            start_of_year = datetime(self.year+1,1,1);
        else:
            start_of_year = datetime(self.year,1,1);
            
        # The end of the file ns41_041226 looks strange at the end as decimal days go negative?!?!?
        # So we take the largest index and start to 
        max_index = np.argmax(doys)

        if(max_index==len(doys)-1):
            return start_of_year + timedelta(days=float(doys[-1]));
        else:
            # Handle the strange end of file decimal days here... 
            doy = doys[max_index];
            new_index = max_index+1;
            while((doy>0) and (new_index<len(doys))):
                doy = doys[new_index];
                new_index += 1;

            if(new_index!=len(doys)):# Something wierd has happened to the decimal days
                return start_of_year + timedelta(days=float(doys[new_index-2]));
            else:
                return start_of_year + timedelta(days=float(doys[-1]));
    
        
    def get_next_datecode(self):
        newdate  = self.get_last_date();
        
        
        # newdate = start_of_year + timedelta(days=float(doy));
        return str(newdate.year).zfill(2) + str(newdate.month).zfill(2) + str(newdate.day).zfill(2);
        
        
    # Return a string with the path to the data on the web
    def folder_path(self):
        # Add check on whether the satellite number is sensible here.
        # The data web address 
        base = 'https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/gps/data/'
        return base + 'ns' + str(self.satellite_number) + '/';

    # To unit test this function 
    # fpath = folder_path(41,07,01,01);
    # print('Folder path is : ', fpath); 
    
    # Return a string with the path to the filename on the web
    def file_path(self):
        # Add input checks here (e.g. is the day less than 31?)  
        # Datafile version number 
        version_number = 'v1.03';# Check whether this varies between the satellite datasets?
        # Create the filename 
        return 'ns' + str(self.satellite_number) + '_' + str(self.year).zfill(2) + str(self.month).zfill(2) + str(self.day).zfill(2) + '_' + version_number + '.ascii';
    # To unit test this function 
    # filepath = file_path(41,07,01,01);
    # print('File path is : ', filepath); 


    # Helper function to find the file size
    def file_size(self,fname):  
            statinfo = os.stat(fname)  
            return statinfo.st_size  


    # Function to download the selected data
    def download(self):

        # Download file
        wget_file = wget.download(self.basepath+self.data_filename);

        # Check if we've downloaded the 404 notice or actual data
        fs = self.file_size(self.data_filename);#This throws an error when wget hasn't finished 
        if(fs<10000):
            print(' - Check the filename as the returned file appears too small');
            print(' - ... deleting file');
            os.remove(self.data_filename);
            return 1;
        else:
            return 0;
    # Unit test : This should download a 404 (not found) file and then delete it (done = 1)
    #done = download_data(41,01,01,01);
    #print(done);
    # Unit test : This should download an actual data file (done = 0)
    #done = download_data(41,07,01,01);
    #print(done);

    def webread(self):
        
        text = urlopen(self.data_filename).read();        
        self.raw_data = text.split('\n');
        
        # Creates a file header to store info about columns for one iteration
        self.json_filename = self.data_filename[-23:-6] + '.json';
        
        # Opens a file for the JSON header
        header = open(self.json_filename, "w+")
        
        # The number of header lines
        header_line_number = 0;        
        
        for i,line in enumerate(text.split('\n')):
            if(i>4):
                if line.startswith('#'):
                    if(i>header_line_number):
                        header_line_number = i;
                    li = line.strip('#');
                    header.write(li+str('\n'));
        header.close();


        # Update self.raw_data as a float array of dimensions Ndata_points by Nparameters
        data_array = np.array([]);
        for i,line in enumerate(self.raw_data):
            if(i>header_line_number):# The data is from here on
                newrow = line.split();
                if(i==header_line_number+1):
                    data_array = np.asarray(newrow);
                else:
                    if(len(newrow)>0):
                        data_array = np.vstack([data_array, newrow]);
 
        self.raw_data = data_array.astype(np.float);


        
        # Open up the json file    
        head = json.load(open(self.json_filename)); #load header file

        # print ' - Reading file : ' + self.data_filename;
        
        variables = {};
        for x, value in head.items():
            i = value['START_COLUMN']
            j = i + value['DIMENSION'][0]
            variables.update({x: self.raw_data[ : , i : j] for column in head});
            
        return variables; 


        # return self.json_filename;



    # Loads a data file using information from the JSON header file 
    def read(self):

        # Open up the json file    
        head = json.load(open(self.json_filename)); #load header file
        f = np.loadtxt(self.data_filename); #use numpy library to create arrays

        print ' - Reading file ' + self.data_filename;
        
        variables = {};
        for x, value in head.items():
            i = value['START_COLUMN']
            j = i + value['DIMENSION'][0]
            variables.update({x: f[ : , i : j] for column in head});
            # print variables;
        
        
        
        return variables;

    # Unit test : Loads the data and prints it
    # data_filename = file_path(41,07,01,01);
    # data = load_data(create_json_header(data_filename),data_filename)
    # print(data);
    
    # Creates a JSON header from the datafile
    def create_json(self):

        # Opens up the datafile
        log = open(self.data_filename, 'r');
        reading = log.readlines()[5:];
    
        # Creates a file header to store info about columns for one iteration
        self.json_filename = self.data_filename[:-6] + '.json';
        header = open(self.json_filename, "w+")
        for line in reading:
            if line.startswith('#'):
                li = line.strip('#')
                header.write(li)
        header.close()
        return self.json_filename;
    # Unit test : Creates a JSON file corresponding to the .ascii file
    # filename = file_path(41,07,01,01);
    # create_json_header(filename);# Should see a JSON file in the directory now


#====================================================================================================
#        
# This class is holds all of the data and functions associated with a single satellite
#
#====================================================================================================
class gps_satellite_data():

    
    def __init__ (self,satellite_number,start_date,end_date):                  

        self.save_disk_space = True;# Deletes the downloaded file once the data has been loaded into this data structure

        # Create a list of datasets     
        self.dataset = [];    
        
        self.start_date = start_date;
        self.end_date = end_date;
        
        self.satellite_number=satellite_number;
        self.filenames = self.get_datafiles_in_date();
        
        if(self.filenames==list()):
            self.empty = True;
            return;
        else:
            self.empty = False;
        
        print; 
        print '====================================';
        print 'Loading data for satellite ', self.satellite_number;

        base = 'https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/gps/data/';
        folder_path = base + 'ns' + str(self.satellite_number) + '/';

        # Load each file 
        for f in self.filenames:

            # Open the particular file 
            full_file_path = folder_path+f;
            print 'File : ' + full_file_path;
            # Append the first datafile into the list 
            self.dataset.append(gps_particle_datafile(full_file_path, self.satellite_number));
    


    def get_datefile_list(self):
        base = 'https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/gps/data/'
        folder_path = base + 'ns' + str(self.satellite_number) + '/';
        urlpath = urlopen(folder_path);
        files = urlpath.read().decode('utf-8');
        
        
        start_file_positions = re.finditer('ns' + str(self.satellite_number), files);
        # end_file_positions = re.finditer('.ascii', files);
        
        # Search for the filenames 
        filelist=[];
        for m in start_file_positions:
            filename = files[m.start():m.start()+23];
            filelist.append(filename);
        return filelist[1::2]; 

    def get_datafiles_in_date(self):
        filelist = self.get_datefile_list();
        cut_list = [];
        
        for f in filelist:
            year = int(f[5:7])+2000;
            month = int(f[7:9]);
            day = int(f[9:11]);
            d = datetime(year,month,day);
            if ( (d >= self.start_date ) and (d <= self.end_date) ):
                cut_list.append(f);
                
        return cut_list;
        

        
    def get_next_datafile(self):

        # Get the next file's datecode 
        new_date_code = self.dataset[-1].get_next_datecode()

        # Read the next file - Load the new data 
        year = int(new_date_code[0:2]);
        month = int(new_date_code[2:4]);
        day = int(new_date_code[4:6]);
        
        self.dataset.append(gps_particle_datafile(self.satellite_number,day,month,year));


        # Check that we were able to load it             
        if(self.dataset[-1].fail_load):
            del self.dataset[-1];# Failed to load so remove it from the list
            return 1;
        else:
            return 0; 
        

#==============================================================================
#         # Get the current size of the dataset list 
#         curr_len = len(self.dataset);
# 
#         # Get the first and last day of the current datafile        
#         first_day = int(math.ceil(self.dataset[-1].summary.min_decimal_day));
#         last_day = int(math.ceil(self.dataset[-1].summary.max_decimal_day));
#             
#         # For the 2004-2005 new year, the decimal day goes negative 
#         # across the end of year boundary ... need to investigate why this is? 
#         if(last_day<0):
#             last_day = 1+int(math.ceil(np.max(self.dataset[-1].data["decimal_day"])));
# 
#         # print first_day;
#         # print last_day;
# 
#         # Convert the decimal day to month and the day
#         epoch = datetime(2000+self.current_year - 1, 12, 31);
#         result = epoch + timedelta(days=last_day);
# 
#         # print result; 
#         # self.dataset[-1].summary.print_summary();
#             
#         # Find the new month and day
#         self.current_day = result.day;
#         self.current_month = result.month;
#         if(first_day>last_day):
#             self.current_year = self.current_year + 1;
#         else:
#             self.current_year = result.year-2000;
#             
#         # Load the new data 
#         self.dataset.append(gps_particle_datafile(self.satellite_number,self.current_day,self.current_month,self.current_year));
# 
#         # Delete the file if this flag is set                
#         if(self.save_disk_space):
#             print 'Saving disk space, deleting file : ', self.dataset[-1].data_filename;
#             os.remove(self.dataset[-1].data_filename);
# 
#         if(curr_len==len(self.dataset)):
#             return 1;
#         else:
#             return 0;
# 
# # UNIT TEST: 
# # This loads the first datafile
# # satellite_data = gps_satellite_data(41,16,12,01);
# # This loads all of the data for the single satellite
# # while(satellite_data.get_next_datafile()==0):
# #    print;
# 
#==============================================================================
# This holds the data associated with an event (e.g. earthquake)
class event():
    def __init__(self,name):
        
        self.name = name;
        
        # Geographic latitude and longitude
        self.data = dict();
        
        self.first_set=True;

    def show_line(self,axes, x_axis_units = 'decimal_year'):
        if(not hasattr(self, 'date')):
            print 'Date property not set.';
            return;
        else:
            if(x_axis_units=='decimal_year'):
                d = self.date;
                dec_year = (float(d.strftime("%j"))-1) / 366 + float(d.strftime("%Y"));
                ylims = axes.get_ylim();
                plt.plot([dec_year,dec_year],ylims,'k--');
            
            if(x_axis_units=='datetime'):
                d = self.date;
                ylims = axes.get_ylim();
                plt.plot([d,d],plt.gca().get_ylim(),'k--');
            
    

    def add_data(self,name,value):
        self.data[name] = value;
        
    def add_date(self,datetime_in):
        self.date = datetime_in;
        self.add_data("year",datetime_in.year);
        
    def add_date(self,mins,hh,dd,mm,yy):
        
        # Save in a date structure 
        if(yy<2000):
            self.date = datetime(day=dd,month=mm,year=2000+yy,hour=hh,minute=mins);
            self.add_data("year",2000+yy);
            startOfThisYear = datetime(year=yy+2000, month=1, day=1);
            startOfNextYear = datetime(year=yy+2000+1, month=1, day=1);
        else:
            self.date = datetime(day=dd,month=mm,year=yy,hour=hh,minute=mins);
            self.add_data("year",yy);
            startOfThisYear = datetime(year=yy, month=1, day=1);
            startOfNextYear = datetime(year=yy+1, month=1, day=1);

        # print self.date;
        year = self.date.year
        

        # secs_so_far = float((self.date - startOfThisYear).seconds);
        # secs_in_year = float((startOfNextYear-startOfThisYear).seconds);
        # print secs_so_far;
        
        decimal_day = float((self.date - startOfThisYear).days) + (float(hh)/24.0) + (float(mins)/(24.0*60.0));
        
        # print decimal_day;
        self.add_data("decimal_day",decimal_day);
        self.add_data("datetime",self.date);
        
        # Save decimal year too 
        # days_in_year = float((startOfNextYear-startOfThisYear).days);       
        # print days_in_year;
        # year_fraction = decimal_day/days_in_year;
        # self.add_data("decimal_year",2000+yy+year_fraction);


        
        
  #      self.add_data("decimal_year",fraction+yy+2000);

        # ayear = (datetime(year=2000+yy+1,month=1,day=1)-datetime(year=2000+yy,month=1,day=1)).total_seconds();
        # td = (datetime(day=dd,month=mm,year=2000+yy,hour=hh,minute=mins) - datetime(year=2000+yy,month=1,day=1)).total_seconds();
        
        # self.add_data("Decimal_",yy);
        # print td/ayear
        # print year_elapsed.year
        

# A basic class to search for earthquakes...
class earthquake_search:
    
    def __init__(self,startdate,enddate,min_magnitude=0,min_lat=-90,max_lat=90,min_lon=-180,max_lon=180):
        
        self.information_text = self.get_earthquake_information(startdate,enddate,min_magnitude,min_lat,max_lat,min_lon,max_lon);
        self.earthquake_information = self.parse_result(self.information_text);
        
    def get_lat_lon(self):
        lats = (self.earthquake_information[:,2]);
        lons = (self.earthquake_information[:,3]);
        return (lats.astype(np.float),lons.astype(np.float));
        
    def get_magnitude(self):
        return self.earthquake_information[:,4];

    def get_datetimes(self):
        
        datetimes = list();
        timestamps = self.earthquake_information;
        for ts in timestamps:
            t = ts[1];
            year = int(t[0:4]);
            month = int(t[5:7]);
            day = int(t[8:10]);
            hour = int(t[11:13]);
            mins = int(t[14:16]);
            secs = int(t[17:19]);
            dt = datetime(year,month,day,hour,mins,secs);
            datetimes.append(dt);
            
        return datetimes;

    # Calculate the L-shell parameters
    def get_L_shells(self, altitude_EM_capture_by_flux_tube_in_km):
        """ 
            Convert geographic latitude to geomagnetic latitude
            See Aleksandrin et al, "High-energy charged particle bursts in the near-Earth space as earthquake precursors", Annales Geophysicae (2003) 21: 597-602
            The altitude of EQ L-shell should be the altitude of the EM wave coupling to the flux tube, estimated to be 300km
        """
        [lat,lon] = self.get_lat_lon();
        mlats=[];
        mlons=[];
        
        L_shells=[];
        for i,la in enumerate(lat):
            [mlat,mlon] = aacgmv2.convert(lat[i], lon[i],altitude_EM_capture_by_flux_tube_in_km);
            r_Re = (6371.0+altitude_EM_capture_by_flux_tube_in_km)/6371.0;
            L_shells.append(r_Re/(math.cos((math.pi/180.0)*mlat)*math.cos((math.pi/180.0)*mlat)));

        return np.asarray(L_shells);

    def get_events(self):
        events = list();
        for ee in self.earthquake_information:
            this_event = event(ee[-1]);

            # Get the date information
            t = ee[1];
            year = int(t[0:4]);
            month = int(t[5:7]);
            day = int(t[8:10]);
            hour = int(t[11:13]);
            mins = int(t[14:16]);
            secs = int(t[17:19]);

            this_event.add_date(mins,hour,day,month,year);

            # Get the lat lon
            this_event.add_data("Geographic_Latitude",float(ee[2]));
            this_event.add_data("Geographic_Longitude",float(ee[3]));

            events.append(this_event);
            
        return events;

    def get_earthquake_information(self,startdate,enddate,min_magnitude=0,min_lat=-90,max_lat=90,min_lon=-90,max_lon=90):

        base_string = 'https://earthquake.usgs.gov/fdsnws/event/1/query?format=text&starttime=';
        mid_string = '&endtime=';
        mag_string = '&minmagnitude=';
        min_lon_string = '&minlongitude=';
        max_lon_string = '&maxlongitude=';
        min_lat_string = '&minlatitude=';
        max_lat_string = '&maxlatitude=';
        
        # Construct string 
        #  e.g. https://earthquake.usgs.gov/fdsnws/event/1/count?starttime=2014-01-01&endtime=2014-01-02
        
        start_time_string = startdate.strftime('%y-%m-%d') + 'T' + str(startdate.hour).zfill(2) + ':' + str(startdate.minute).zfill(2) + ':' + str(startdate.second).zfill(2);
        end_time_string = enddate.strftime('%y-%m-%d') + 'T' + str(enddate.hour).zfill(2) + ':' + str(enddate.minute).zfill(2) + ':' + str(enddate.second).zfill(2);
        
        url_string = base_string + start_time_string  + mid_string + end_time_string + mag_string + str(min_magnitude) + min_lon_string + str(min_lon) + max_lon_string + str(max_lon) + min_lat_string + str(min_lat) + max_lat_string + str(max_lat); 
        # print url_string;
        
        self.information_text = urlopen(url_string).read();
        self.earthquake_information = self.parse_result(self.information_text)
        
        return self.information_text;
        
    
    # This function was copied from  http://qingkaikong.blogspot.co.uk/2016/02/query-usgs-catalog-online-and-plot.html
    def parse_result(self,inputText):
        '''
        Function to parse the requested earthquake events data from USGS, and save it into numpy array
        '''
        event_id = []
        origin_time = []
        evla = []
        evlo = []
        evdp = []
        mag = []
        mag_type = []
        EventLocationName  = []
        for i, item in enumerate(inputText.split('\n')[0:-1]):
            if i < 1:
                # skip the header
                continue
                
            try:
                splited = item.split('|')
                event_id.append(splited[0])  
                origin_time.append(splited[1])
                evla.append(splited[2])
                evlo.append(splited[3])
                evdp.append(splited[4])
                mag.append(splited[10])
                mag_type.append(splited[9])
                EventLocationName.append(splited[-1])
            except:
                # just in case there are some wrong data in the catlog
                print item
                print 'Skip wrong data or there is something wrong' 
        
        return np.c_[event_id, origin_time, evla, evlo, mag, mag_type, EventLocationName]












class satellite_data_plot:

    
    
    
    def __init__(self):
        self.x_axis_data='decimal_days';
        self.marker_colour = 'b';
        self.marker_size = 3;        
        
        
    def daterange(start_date, end_date, step_in_days):
        days_range = np.arange(0,(end_date - start_date).days, step_in_days)
        for n in days_range:
            yield start_date + timedelta(days=n);

    def get_counts(self,
                   output_data,
                   satellite_label = 'ns41',
                   data_channel = 0,
                   keyvarname = 'rate_proton_measured'):
               
        meas_period = np.asarray(output_data[satellite_label]['collection_interval']);
        rate = np.asarray(output_data[satellite_label][keyvarname])[:,data_channel];        
        
        return np.multiply(meas_period,rate);
        
    def get_signal_to_background_ratio(self,
                   output_data,
                   satellite_label = 'ns41',
                   data_channel = 0,
                   keyvarname = 'rate_proton_measured'):

        signal = np.asarray(output_data[satellite_label][keyvarname])[:,data_channel];

        if('proton' in  keyvarname):
            bkgname = 'proton_background';
        if('electron' in  keyvarname):
            bkgname = 'electron_background';            
        background = np.asarray(output_data[satellite_label][bkgname])[:,data_channel];        
        
        return np.divide(signal,background);
    
    def get_signal_to_stddev_ratio(self,
                   output_data,
                   satellite_label = 'ns41',
                   data_channel = 0,
                   keyvarname = 'proton'):

        signal = np.asarray(output_data[satellite_label][keyvarname])[:,data_channel];

        if('proton' in  keyvarname):
            bkgname = 'proton_background';
        if('electron' in  keyvarname):
            bkgname = 'electron_background';            
        stddev_bkg = np.asarray(math.sqrt(output_data[satellite_label][bkgname][:,data_channel]));        
        
        return np.divide(signal,stddev_bkg);
    

    def get_li_ma_significance(self,
                               output_data,
                               on_startdate,
                               on_enddate,
                               off_startdate,
                               off_enddate,
                               satellite_label = 'ns41',
                               keyvarname = 'rate_proton_measured', 
                               data_channel = 0):
        
        # First get the number of particles 
        counts = self.get_counts(output_data, satellite_label, data_channel,keyvarname);
        
        # Calculate the on and off counts 
        on_date_cuts = self.get_date_cuts(counts,on_startdate,on_enddate);
        on_counts = sum(counts[on_date_cuts]);

        off_date_cuts = self.get_date_cuts(counts,off_startdate,off_enddate);
        off_counts = sum(counts[off_date_cuts]);
        
        # Find the ratio of the time periods
        alpha = (on_enddate - on_startdate).total_seconds()/(off_enddate - off_startdate).total_seconds();

        # Return the significance 
        num1 = (1+alpha)*on_counts;
        den1 = alpha*(on_counts + off_counts);
        num2 = (1+alpha)*off_counts;
        den2 = on_counts + off_counts;
        S = math.sqrt(2)*math.sqrt( on_counts * math.log(num1/den1) + off_counts*math.log(num2/den2) );        

        return S;                               
    

    def plot_li_ma(self, 
                   fig_no,
                   axes,
                   output_data, 
                   satellite_label='ns41',
                   keyvarname = 'rate_proton_measured',     
                   data_channel = 0,
                   startdate=datetime(1990,1,1),
                   enddate=datetime(2017,1,1), 
                   on_off_period_days = 1,
                   Lmin = 0, Lmax = 100):

        # Loop over the range of days, with         
        for this_date in self.daterange(startdate,enddate,on_off_period_days):
            print this_date;
        
    

        #==============================================================================
        #   Simply plots the data given the date and L-shell cuts 
        #==============================================================================
    def plot_channels(self,fig_no,axes,output_data, 
                     satellite_label='ns41',
                     keyvarname = 'rate_proton_measured',     
                     data_channel = 0,
                     startdate=datetime(1990,1,1),enddate=datetime(2017,1,1), 
                     Lmin = 0, Lmax = 100):

        fig = plt.figure(fig_no);
        colour = ['red', 'blue', 'green', 'magenta', 'navy', 'aqua', 'orange', 'yellow' , 'lime', 'peru', 'slategrey', 'teal', 'maroon', 'mediumspringgreen', 'red']

        # The x-axis data
        x =  np.asarray(output_data[satellite_label]['datetime']);                    
        # x = self.daterange(startdate,enddate,0.25);
        
        # The signal data            
        y = np.asarray(output_data[satellite_label][keyvarname])[:,data_channel];

        # Get the date cuts 
        this_data = np.squeeze(np.asarray(output_data[satellite_label]["datetime"]));
        date_cuts = self.get_date_cuts(this_data,startdate,enddate);

        # And the L shell cuts 
        Lshell_data = np.squeeze(np.asarray(output_data[satellite_label]["L_shell"]));
        Lshell_cuts = self.get_Lshell_cuts(Lshell_data,Lmin,Lmax);

        # Make sure none of the data has been dropped 
        dropped = np.squeeze(np.asarray(output_data[satellite_label]["dropped_data"]));
        dropped_cuts = self.get_dropped_data_cuts(dropped);
        
        data_cuts = (date_cuts & Lshell_cuts & dropped_cuts );

        plt.plot(x[data_cuts],y[data_cuts],'o-',color=self.marker_colour);#,s=self.marker_size);
 
        # axes.legend(loc = 'upper left');
        plt.xlabel('Date');
        plt.xticks(rotation=50);
        plt.ylabel((keyvarname.replace("_", " ")).title());
        plt.xlim(startdate,enddate);       
    
    
    def get_date_cuts(self,datetimes,startdate,enddate):
        date_cuts=[];
        for idate in datetimes:
            date_cuts.append( (idate>startdate) and (idate<enddate) );
        
        return np.asarray(date_cuts);

    
    def get_Lshell_cuts(self,Lshells,Lmin,Lmax):
        Lshell_cuts=[];
        for iL in Lshells:
            Lshell_cuts.append( (iL>Lmin) and (iL<Lmax) );
        return np.asarray(Lshell_cuts);
    
    
    def get_dropped_data_cuts(self,dropped):
        dropped_cuts=[];
        for d in dropped:
            dropped_cuts.append( d == 0 );
        return np.asarray(dropped_cuts);
    
    









    
class satellite_map_plot:
    def __init__(self):

        # miller projection
        self.map_plot = Basemap(projection='mill',llcrnrlon=60,llcrnrlat=-60, urcrnrlon=150, urcrnrlat=60)
       
        # plot coastlines, draw label meridians and parallels.
        self.map_plot.drawcoastlines()
        self.map_plot.drawmapboundary(fill_color='aqua')
        self.map_plot.drawcountries(color='0.4')
        self.map_plot.fillcontinents(color='coral',lake_color='aqua')
        self.map_plot.drawmapboundary(fill_color='0.8')
        self.map_plot.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0]);

        self.map_plot.drawmeridians(np.arange(self.map_plot.lonmin,self.map_plot.lonmax+30,60),labels=[0,0,0,1]);
        
        # Formatting parameters 
        #  - For the events 
        self.marker_scale=1.0;
        self.marker_colour='b';
        self.inner_colour='r';
        self.outer_colour='g';
        self.data_keyname='';
        
        self.label_lat = -87.5;
        self.label_lon = -80;

    def plot_event_on_map(self,fig_no,events):
    
        plt.figure(fig_no);

        # fill continents 'coral' (with zorder=0), color wet areas 'aqua'
        # shade the night areas, with alpha transparency so the
        # map shows through. Use current time in UTC.
        date = datetime.utcnow()
        for e in events:
            # Plot the Earthquake event on the map
            lon = np.squeeze(np.asarray(e.data["Geographic_Longitude"]));
            lat = np.squeeze(np.asarray(e.data["Geographic_Latitude"]));
            # print lat,lon
            x,y = self.map_plot(lon,lat);
            self.map_plot.scatter(x,y,s=300*self.marker_scale,marker='o',color=self.outer_colour,alpha=0.50,zorder=5);
            self.map_plot.scatter(x,y,s=50*self.marker_scale,marker='o',color=self.inner_colour,alpha=0.5,zorder=5);
            self.map_plot.scatter(x,y,s=1000*self.marker_scale,marker='x',color='k',alpha=1,zorder=5);

    def signal_heatmap(self,output_data,fig_no, axes, keyvarname = 'rate_proton_measured',
                       data_channel = 0,
                       satellite_labels=['ns41'],startdate=datetime(1990,1,1),enddate=datetime(2017,1,1), Lmin = 0, Lmax = 100):

        fig = plt.figure(fig_no);
        
        # Get the output data by satellite
        # output_data = self.get_selected_data_by_satellite();

        xs=[];
        ys=[];
        zs=[];
        for sat in satellite_labels:
            # Get the lat and long 
            sat_lat,sat_lon = self.get_lat_long(output_data, sat);
            xs1,ys1 = self.map_plot(sat_lon,sat_lat);
            
            # Get the date cuts 
            this_data = np.squeeze(np.asarray(output_data[sat]["datetime"]));
            date_cuts = self.get_date_cuts(this_data,startdate,enddate);

            # And the L shell cuts 
            Lshell_data = np.squeeze(np.asarray(output_data[sat]["L_shell"]));
            Lshell_cuts = self.get_Lshell_cuts(Lshell_data,Lmin,Lmax);

            # Make sure none of the data has been dropped 
            dropped = np.squeeze(np.asarray(output_data[sat]["dropped_data"]));
            dropped_cuts = self.get_dropped_data_cuts(dropped);
        
            data_cuts = (date_cuts & Lshell_cuts & dropped_cuts );

            # Get the colour data 
            signal_data = np.squeeze(np.asarray(output_data[sat][keyvarname]));
            signal_data = signal_data[:,int(data_channel)]; 

            xs1 = np.ravel(xs1[data_cuts]);
            ys1 = np.ravel(ys1[data_cuts]);
            zs1 = np.ravel(signal_data[data_cuts]);

            for x in xs1:
                xs.append(x);
            for y in ys1:
                ys.append(y);
            for z in zs1:
                zs.append(z);

        xs = np.asarray(xs);
        ys = np.asarray(ys);
        zs = np.asarray(zs);
        
        xmin = xs.min()
        xmax = xs.max()
        ymin = ys.min()
        ymax = ys.max()

        gridsize=60;
        
        hb = plt.hexbin(xs, ys, C=zs, gridsize=gridsize, marginals=True, cmap=plt.cm.PuRd,     zorder=100);

        cb = fig.colorbar(hb, ax=axes)
        cb.set_label('Satellite Visits')


#==============================================================================
# # http://bagrow.com/dsv/heatmap_basemap.html
#     def signal_heatmap(self,output_data):
#                 
#         
#                 # ######################################################################
#                 # bin the epicenters (adapted from 
#                 # http://stackoverflow.com/questions/11507575/basemap-and-density-plots)
#         
#             # compute appropriate bins to chop up the data:
#                 db = 1 # bin padding
#         lon_bins = np.linspace(min(lons)-db, max(lons)+db, 10+1) # 10 bins
#         lat_bins = np.linspace(min(lats)-db, max(lats)+db, 13+1) # 13 bins
#             
#         density, _, _ = np.histogram2d(lats, lons, [lat_bins, lon_bins])
#         
#         # Turn the lon/lat of the bins into 2 dimensional arrays ready
#         # for conversion into projected coordinates
#         lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
#         
#         # convert the bin mesh to map coordinates:
#         xs, ys = m(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh
#         # ######################################################################
# 
#==============================================================================
    def show_satellite_label(self,fig_no, output_data, satellite_label='ns41'):

        
        if(self.label_lon < -90):
            self.label_lon+=360;
        #        sat_lat,sat_lon = self.get_lat_long(output_data, satellite_label);
        xs,ys = self.map_plot(self.label_lon,self.label_lat);
        
        plt.text(xs, ys, satellite_label, fontsize=8, fontweight='bold',bbox=dict(facecolor='w', edgecolor=self.marker_colour, linewidth=2, alpha=0.85), zorder=100);
        self.label_lat += 0;
        self.label_lon += 30;
            
            
    def get_normalized_data(self, output_data, satellite_label, keyvarname):            

        data = np.squeeze(np.asarray(output_data[satellite_label][keyvarname][:]))[:,0];
        data_scaled = np.zeros(data.shape);
        for index, d in enumerate(data):
            if(d!=-1):
                data_scaled[index] = (d - min(data)) / -np.ptp(data);#(max(data) - min(data));
            else:
                data_scaled[index] = 0;
        
        return data_scaled;
                

    def get_date_cuts(self,datetimes,startdate,enddate):
        date_cuts=[];
        for idate in datetimes:
            date_cuts.append( (idate>startdate) and (idate<enddate) );
        
        return np.asarray(date_cuts);

    
    def get_Lshell_cuts(self,Lshells,Lmin,Lmax):
        Lshell_cuts=[];
        for iL in Lshells:
            Lshell_cuts.append( (iL>Lmin) and (iL<Lmax) );
        return np.asarray(Lshell_cuts);
    
    
    def get_dropped_data_cuts(self,dropped):
        dropped_cuts=[];
        for d in dropped:
            dropped_cuts.append( d == 0 );
        return np.asarray(dropped_cuts);
    
    
    def get_lat_long(self,output_data, satellite_label):
        sat_lon = np.asarray(output_data[satellite_label]["Geographic_Longitude"]);
        sat_lon[sat_lon < -90]+=360;
        sat_lat = np.asarray(output_data[satellite_label]["Geographic_Latitude"]);
        return sat_lat,sat_lon;
    
    # Plot the projection of the satellite
    def plot_satellite_trajectory_colour_by_data(self,output_data,fig_no, keyvarname = 'rate_proton_measured',                              
                              satellite_label='ns41',startdate=datetime(1990,1,1),enddate=datetime(2017,1,1), Lmin = 0, Lmax = 100):

        plt.figure(fig_no);
        
        # Get the output data by satellite
        # output_data = self.get_selected_data_by_satellite();

        # Get the lat and long 
        sat_lat,sat_lon = self.get_lat_long(output_data, satellite_label);
        xs,ys = self.map_plot(sat_lon,sat_lat);

        # Get the date cuts 
        this_data = np.squeeze(np.asarray(output_data[satellite_label]["datetime"]));
        date_cuts = self.get_date_cuts(this_data,startdate,enddate);

        # And the L shell cuts 
        Lshell_data = np.squeeze(np.asarray(output_data[satellite_label]["L_shell"]));
        Lshell_cuts = self.get_Lshell_cuts(Lshell_data,Lmin,Lmax);

        # Make sure none of the data has been dropped 
        dropped = np.squeeze(np.asarray(output_data[satellite_label]["dropped_data"]));
        dropped_cuts = self.get_dropped_data_cuts(dropped);
        
        data_cuts = (date_cuts & Lshell_cuts & dropped_cuts );

        # Colour by the signal data            
        data_scaled = self.get_normalized_data(output_data, satellite_label, keyvarname);            
        self.map_plot.scatter(xs[data_cuts],ys[data_cuts],
                                 s=25*self.marker_scale,c=data_scaled[data_cuts],marker='o',zorder=100,cmap='jet',alpha=1);

    # Plot the projection of the satellite
    def plot_satellite_trajectory(self,output_data,fig_no,                               
                              satellite_label='ns41',startdate=datetime(1990,1,1),enddate=datetime(2017,1,1), Lmin = 0, Lmax = 100):

        plt.figure(fig_no);
        
        # Get the output data 
        # output_data = self.get_selected_data('date');
        # output_data = self.get_selected_data_by_satellite();

        label_plotted = False;
        colour = ['red', 'blue', 'green', 'magenta', 'navy', 'aqua', 'orange', 'yellow' , 'lime', 'peru', 'slategrey', 'teal', 'maroon', 'mediumspringgreen', 'red']

        # Get the lat and long 
        sat_lat,sat_lon = self.get_lat_long(output_data, satellite_label);
        xs,ys = self.map_plot(sat_lon,sat_lat);

        # Get the date cuts 
        this_date = np.squeeze(np.asarray(output_data[satellite_label]["datetime"]));
        date_cuts = self.get_date_cuts(this_date,startdate,enddate);

        # And the L shell cuts 
        Lshell_data = np.squeeze(np.asarray(output_data[satellite_label]["L_shell"]));
        Lshell_cuts = self.get_Lshell_cuts(Lshell_data,Lmin,Lmax);
        
        # Make sure none of the data has been dropped 
        dropped = np.squeeze(np.asarray(output_data[satellite_label]["dropped_data"]));
        dropped_cuts = self.get_dropped_data_cuts(dropped);
        
        data_cuts = (date_cuts & Lshell_cuts & dropped_cuts )
        
        # Plot the data 
        self.map_plot.scatter(xs[data_cuts],ys[data_cuts],s=25*self.marker_scale,marker='o',zorder=100,c=self.marker_colour,alpha=0.5);


#==============================================================================
#         for fn in output_data.keys():
#             if(satellite_label==fn):
#                 
#                 if(data_keyname==''):
#                             
#                     sat_lon = np.squeeze(np.asarray(output_data[fn]["Geographic_Longitude"])).flatten();
#                     sat_lon[sat_lon < -90]+=360;
#                     sat_lat = np.squeeze(np.asarray(output_data[fn]["Geographic_Latitude"])).flatten();
# 
#                                         
#                     # Date cut 
#                     this_date = np.squeeze(np.asarray(output_data[fn]["datetime"]));
# 
#                     xs,ys = map(sat_lon,sat_lat);
#                     date_cuts=[];
#                     for i,idate in enumerate(this_date):
#                         date_cuts.append( (idate>startdate) and (idate<enddate) );
#                         
#                     #map.plot(xs,ys,marker='o',color='g',markersize=15);
#                     #map.plot(xs,ys,marker='o',color= colour[i],markersize=3,linewidth=0);
# 
#                     map.scatter(xs[date_cuts],ys[date_cuts],s=16*marker_scale,marker='o',color=marker_color,zorder=10);
#                     
#                         
#==============================================================================
#==============================================================================
#                 elif(data_keyname=='proton_signal_to_background'):
#                     
#                     # Colour by the signal data            
#                     sig_data = np.squeeze(np.asarray(output_data[fn]['rate_proton_measured'][:]))[:,0];
#                     bkg_data = np.squeeze(np.asarray(output_data[fn]['proton_background'][:]))[:,0];
#                     data = sig_data/bkg_data
#                     
#                     data_scaled = np.zeros(data.shape);
#                     for index, d in enumerate(data):
#                         data_scaled[index] = (d - min(data)) / -np.ptp(data);#(max(data) - min(data));
# 
#                     sat_lon = np.asarray(output_data[fn]["Geographic_Longitude"]);
#                     sat_lon[sat_lon < -90]+=360;
#                     sat_lat = np.asarray(output_data[fn]["Geographic_Latitude"]);
#                     xs,ys = map(sat_lon,sat_lat);
# 
#                     # Date cut 
#                     this_date = np.squeeze(np.asarray(output_data[fn]["datetime"]));
# 
#                     xs,ys = map(sat_lon,sat_lat);
#                     date_cuts=[];
#                     for i,idate in enumerate(this_date):
#                         date_cuts.append( (idate>startdate) and (idate<enddate) );
# 
#                     map.scatter(xs[date_cuts],ys[date_cuts],s=25,c=data_scaled[date_cuts],marker='o',zorder=100,cmap='jet',alpha=1);
# 
#                     # Indicate the time 
#                     n = int((enddate-startdate).total_seconds()/3600.0)*3;
#                     time_scales = np.zeros(int(len(date_cuts)),dtype=int);
#                     time_scales[0::n] = 1;
#                     map.scatter(xs[date_cuts],ys[date_cuts],s=4*25*time_scales,c='k',marker='s',zorder=50,cmap='jet',alpha=1);
# 
#==============================================================================
                    #if(not label_plotted):
                    #    plt.text(xs[0], ys[0], fn[0:4], fontsize=16, fontweight='bold',bbox=dict(facecolor='w', alpha=0.5), zorder=100);
                    #    label_plotted=True;




# This holds all of the data that searches over all satellites
class meta_search:
    def __init__(self):
        self.satellites = list();
        #self.satellites.append(41);
        #self.satellites.append(48);
	self.satellites.append(61)
	#self.satellites.append(54)
        #for i in np.arange(53,74): 
            #self.satellites.append(i);
        
        self.satellites_data = list();
        self.searches = list();

        
    def plot_all_channels(self,axes,
                  satellite_numbers = ['ns41'], 
                  signal_label = 'rate_proton_measured',
                  background_label = 'proton_background', 
                  x_axis_data='decimal_year'):

        # Some default colours 
        colour = ['rosybrown', 'blue', 'red', 'green', 'magenta', 'navy', 'aqua', 'orange', 'yellow' , 'lime', 'peru', 'slategrey', 'teal', 'maroon', 'mediumspringgreen']

        # Get the output data 
        output_data = self.get_selected_data('date');
        
        
        # Loop over the keys and match 
        data_index = 0;
        for sat in satellite_numbers:
            # print 'Plotting satellite ', sat
            for key in output_data.keys():
                if(sat in key):
                    # print ' - key : ', key;
            
                    # The x-axis data
                    x =  output_data[key][x_axis_data][:];
                    
                    # The signal data            
                    tmp = output_data[key][signal_label][:];

                    k = len(tmp[0])
                    p_eng = np.array_split(tmp, range(k), axis = 1)
                
                    # Convert decimal years to day of year 
                    for j in range(k):
                        axes.scatter(x, p_eng[j+1], color = colour[j])

        # axes.legend(loc = 'upper left');
        plt.xlabel('Decimal Year');
        plt.ylabel((signal_label.replace("_", " ")).title());

    # Searches for a single event over all satellites between the dates
    def apply_search(self,this_event, tol=0):
        for sats_data in self.satellites_data:
            self.searches.append(search(this_event,sats_data,tol));


    def get_signal_to_stddev_ratio(self,
                   output_data,
                   satellite_label = 'ns41',
                   data_channel = 0,
                   keyvarname = 'proton'):

        if('proton' in  keyvarname):
            signame = 'rate_proton_measured'
            bkgname = 'proton_background';
        if('electron' in  keyvarname):
            signame = 'rate_electron_measured'
            bkgname = 'electron_background';            

        signal = (np.asarray(output_data[satellite_label][signame])[:,:,data_channel]).flatten();
        stddev_bkg = ((np.asarray( output_data[satellite_label][bkgname] )[:,:,data_channel]).flatten());
        
        return np.divide( signal, stddev_bkg);

    def get_burst_events(self, output_data):
        
        # output_data = self.get_selected_data_by_satellite();
        # print output_data.keys();
        
        
        burst_events=[];
        
        for sat in output_data:            
                
            signame = 'rate_electron_measured';
            bkgname = 'electron_background';
            
            if(len(output_data[sat].keys())>0):
                signal = np.asarray(output_data[sat][signame])[:,2];
                bkg = np.asarray( output_data[sat][bkgname])[:,2];

		avg = np.mean(signal)
		stddev = np.std(signal)
		sig_dif = np.subtract(signal, avg)
                ratio = np.divide( sig_dif, stddev);
                
                # 
                burst_indices = ratio>4; 
                
                # Get the dates of these events 
                burst_radii = np.asarray(output_data[sat]['Rad_Re'])[burst_indices];
                burst_dates = np.asarray(output_data[sat]['datetime'])[burst_indices];
                burst_latitudes = np.asarray(output_data[sat]['Geographic_Latitude'])[burst_indices];
                burst_longitudes = np.asarray(output_data[sat]['Geographic_Longitude'])[burst_indices];
                burst_Ls = np.asarray(output_data[sat]['L_shell'])[burst_indices];
                
                for i,bd in enumerate(burst_dates):
                    this_event = event('Burst #' + str(i));

                    # Get the date information
                    this_event.add_date(bd.minute,bd.hour,bd.day,bd.month,bd.year);

                    # Get the lat lon
                    this_event.add_data("Geographic_Latitude",float(burst_latitudes[i]));
                    this_event.add_data("Geographic_Longitude",float(burst_longitudes[i]));
                    this_event.add_data("L_shell",float(burst_Ls[i]));
                    this_event.add_data("Rad_Re",float(burst_radii[i]));
                    this_event.add_data("Satellite_Number",float(sat));
                    
                    burst_events.append(this_event);
                    
        return burst_events;


                
                
            
            
            # stddev_ratio = self.get_signal_to_stddev_ratio(output_data, satellite_label=sat, data_channel=3,keyvarname='proton');

#            bursts = stddev_ratio > 5; 
            
  #          for b in bursts:
 #               


   #         print str(sat) + ' ' + str(sum(stddev_ratio>5))+ ' ' + str(len(stddev_ratio)) + ' ' + str(stddev_ratio.shape);
        
    def get_all_data_by_satellite(self):
        
        out_data = dict();
        
        # Loop over all satellites' data
        for sat in self.satellites:
            
            #Initialize data if this is the first time
            if(not (sat in out_data)):
                out_data[sat] = dict();

            # Loop over the data variable names and save to the output data
            for this_dataset_collection in self.satellites_data:# All files
                
                # Check the same satellite label 
                if(sat==this_dataset_collection.satellite_number):
                    
                    # Loop over all data files with this satellite number 
                    for this_data in this_dataset_collection.dataset:
                    
                        # Loop over the list of variables within this file 
                        for var in this_data.data:
                            if(not (var in out_data[sat])):
                                out_data[sat][var] = list();

                            # Save out to the structure 
                            for d in this_data.data[var]:
                                out_data[sat][var].append(d);
                       
                
                    # print sat, this_dataset_collection.satellite_number, (np.asarray(out_data[sat]['decimal_day'])).shape
                    for idd,dd in enumerate(out_data[sat]['decimal_day']):
                    
                        if(not ('datetime' in out_data[sat])):
                            out_data[sat]['datetime'] = list();

                    
                        # Make the datetime object
                        data_year = out_data[sat]['year'][idd];
                        epoch = datetime(int(data_year[0]) - 1, 12, 31);
                        result = epoch + timedelta(days=float(dd));
                        out_data[sat]['datetime'].append(result);
                    
        return out_data; 
        
            
    def get_selected_data_by_satellite(self):
        out_data = dict();
        for s in self.searches:
            for fn in s.selected_data:
                satellite_number = fn[0:4];                         
                
                if(not (satellite_number in out_data)):
                    out_data[satellite_number] = dict();
                
                for var in s.selected_data[fn]:

                    if(not (var in out_data[satellite_number])):
                        out_data[satellite_number][var] = list();
                        
                    for sdata in s.selected_data[fn][var]:
                        out_data[satellite_number][var].append(sdata);

                # Write out the datetime to the output data too (for convenience)
                if(not ('datetime' in out_data[satellite_number])):
                    out_data[satellite_number]['datetime'] = list();
                    
                for idd,dd in enumerate(s.selected_data[fn]['decimal_day']):
                    
                    # Find the datetime object
                    data_year = s.selected_data[fn]['year'][idd];
                    epoch = datetime(data_year - 1, 12, 31);
                    result = epoch + timedelta(days=float(dd));
                    out_data[satellite_number]['datetime'].append(result);
        
        return out_data;

    # Takes the indices and adds extras within a wider time window
    def extend_time_window(self,time_window_in_days):
        for s in self.searches:
            s.get_indices_in_time_window(time_window_in_days);
            
    def clear_searches(self):
        self.searches=list();

    # Loads the data between the dates into the current search         
    def load_data(self,start_date, end_date):
        
        for this_sat in self.satellites:
            # Extract the day,month,year from the start and end dates
            new_data = gps_satellite_data(this_sat,start_date,end_date);
            if(not new_data.empty):
                self.satellites_data.append(new_data);

            # Check that we were able to load it             
            #if(self.satellites_data[-1].dataset[-1].fail_load):
             #   del self.satellites_data[-1];# Failed to load so remove it from the list


#==============================================================================
#     def show_timescales_on_mapself(fig_no,satellite_label='ns41', data_keyname='', 
#                                    marker_scale=1.0, marker_color='b', new_plot=False, startdate=datetime(1990,1,1),enddate=datetime(2017,1,1) ):
#         
#         plt.figure(fig_no);
# 
#         if((not hasattr(self, 'map_plot')) or (new_plot)):
#             map = Basemap(projection='mill',lon_0=90)
#             self.map_plot = map 
#         else:
#             map = self.map_plot;
#         
#         # Get the lat and lon, and datetime objects
#         sat_lon = np.squeeze(np.asarray(output_data[fn]["Geographic_Longitude"])).flatten();
#         sat_lon[sat_lon < -90]+=360;
#         sat_lat = np.squeeze(np.asarray(output_data[fn]["Geographic_Latitude"])).flatten();
#                                         
#         # Date cut 
#         this_date = np.squeeze(np.asarray(output_data[fn]["datetime"]));
# 
#         xs,ys = map(sat_lon,sat_lat);
#         
#         # Indicate the time 
#         n = int((enddate-startdate).total_seconds()/3600.0)*3;
#         time_scales = np.zeros(int(len(date_cuts)),dtype=int);
#         time_scales[0::n] = 1;
#         
#         # Get the local perpendicular at each point 
#         xs_diff = [xs[i]-xs[i-1] for i in range(1,len(xs))];
#         ys_diff = [ys[i]-ys[i-1] for i in range(1,len(ys))];
#    
#         # Eqn. of perpendicular is: yp = perp_slope*xp + cp
#         perp_slope = np.divide(-ys_diff,xs_diff);
#         
#         
#         
#         # Draw a line over the data every so often 
#         # x1 = x0 + 
#         # y1 = y0 + 
#         len_line = 10;
#         for i,g in enumerate(m):               
#         if(time_scales[i]==1): 
# 
#             x1 = xs[i]*perp_slope
#             x2 = xs[i] + t*
# 
#             map.plot();  
#         
#         # map.scatter(xs[date_cuts],ys[date_cuts],s=4*25*time_scales,c='k',marker='s',zorder=50,cmap='jet',alpha=1);
# 
#     
#  
#             
#                 
#                 
#==============================================================================
    def add_keyname_search(self, keyvarnames, tol=0):
        for s in self.searches:
            s.compare(keyvarnames, tol);

    # Returns the data selected by the variable name 
    def get_selected_data(self, key_varname):
        output_data = dict();
        
        if(self.searches==list()):
            print "Error! meta-search::get_selected_data() - No Searches Performed"
            return 1;

        # Loop over all applied searches, each search is the same satellite but over time 
        for isearch,s in enumerate(self.searches):         
            # print isearch, len(self.searches), len(self.satellites_data);
            
            # print s.indices.keys();
            # print s.key_filenames; 
            
            # Get the filenames
            # for fn in s.key_filenames:
            for fn in s.indices.keys():
                
                # Check if the indices set is empty 
                # print "Filename : ", fn, " Empty : ", s.indices[fn][key_varname]!=set(); 
                if(s.indices[fn][key_varname]!=set()):

                    for dataset in self.satellites_data[isearch].dataset:
                        if(dataset.data_filename[0:-12]==fn):
                            output_data[fn] = dict();
                            for ikey in dataset.data.keys():
                                output_data[fn][ikey] = list();
                                for index in s.indices[fn][key_varname]:
                                    output_data[fn][ikey].append(dataset.data[ikey][index]);

        self.selected_data = output_data;                    
        return output_data;        
        
        
    
#==============================================================================
# d = datetime(4,12,26);
# print d.date
# 
#                 
#     
# ms = meta_search();
# 
# e = event("Boxing Day Earthquake");
# e.add_date(mins=59,hh=00,dd=26,mm=12,yy=04);
# e.add_data("Geographic_Latitude",3.24);
# e.add_data("Geographic_Longitude",95.8);
# 
# start_date = datetime(04,12,19);
# end_date = datetime(05,1,1);
# 
# ms.load_data(start_date,end_date);
# ms.apply_search(e);
# ms.extend_time_window(0.5);
# 
# output_data = ms.get_selected_data("decimal_day");
# 
# for fn in output_data.keys():
#     print;
#     print;
#     print fn;
#     for ikey in output_data[fn].keys():
#         print ikey;
#         print output_data[fn][ikey][0];
#         print output_data[fn][ikey][-1];
# 
# 
#==============================================================================
#==============================================================================
# satellites_data = list();
# 
# searches = list();
# for this_sat in sats:
#     
# 
#     # satellite_data.save_disk_space=True;
#     for i in np.arange(0,5):
#         satellites_data[-1].get_next_datafile();
# 
#     # e.add_data("Rad_Re",130);
# 
# 
#     searches.append(search(e,satellites_data[-1]));
#     print "==================================================================="
#     searches[-1].compare_with_zerocrossing(["Geographic_Latitude"]);
#     print "==================================================================="
#     searches[-1].get_indices_in_time_window(0.5);
#     searches[-1].intersect_indices_set();
#     print "==================================================================="
#     searches[-1].print_indices();
#     searches[-1].remove_empty_sets();
#     searches[-1].print_indices();
#     
#==============================================================================








        
        









# Defines a data search 
class search():
    def __init__(self, event, satellite_data,tol=0):
        
        # Create the indices of the satellites
        # self.satellites = np.arange(53,74);
        # self.satellites = np.insert(satellites,0,41);
        # self.satellites = np.insert(satellites,1,48);
        self.satellite_data = satellite_data;
        self.event = event;
        self.indices = dict();# dictionary of sets of indices of the dataset, keyed by filename
        self.selected_data = dict();# dictionary of sets of data
        self.key_varnames = set();
        self.key_filenames = set();
        
        # Initialise indices array
        for sat_data in satellite_data.dataset[:]:
            # Create new set in the indices dictionary
            key_filename = sat_data.data_filename[:-12];
            self.key_filenames.add(key_filename);
            self.indices[key_filename] = dict();
            self.selected_data[key_filename] = dict();

        self.compare_dates(tol);

    def remove_empty_sets(self):
        # Print indices         
        for key_filename in self.key_filenames:
            for key_varname in self.key_varnames:
                if(self.indices[key_filename][key_varname]==set()):
                    print "Deleting " + key_filename + " " + key_varname;
                    del self.indices[key_filename][key_varname];
        
        
    def get_indices_in_time_window(self, time_window_in_days):
        
        # For each search variable 
        for sat_data in self.satellite_data.dataset[:]:
 
            # Create new set in the indices dictionary
            key_filename = sat_data.data_filename[:-12];
            
            new_indices = set();
            for index in self.indices[key_filename]["decimal_day"]:

                decimal_day_event = self.event.data["decimal_day"];
                decimal_day_data = sat_data.data["decimal_day"][index];
                    
                upper_index=0;
                while(((upper_index+index)>0) and (math.fabs(decimal_day_data - decimal_day_event)<time_window_in_days)):
                    # print index, upper_index, decimal_day_data, decimal_day_event;
                    if(not ((index+upper_index) in self.indices[key_filename]["decimal_day"])):
                        new_indices.add(index+upper_index);                    
                    
                    upper_index-=1;
                    decimal_day_data = sat_data.data["decimal_day"][index+upper_index];


                upper_index=0;
                decimal_day_data = sat_data.data["decimal_day"][index];
                while (math.fabs(decimal_day_data - decimal_day_event)<time_window_in_days):
                    # print index,index+upper_index, decimal_day_data, decimal_day_event;
                    
                    if(not ((index+upper_index) in self.indices[key_filename]["decimal_day"])):
                        new_indices.add(index+upper_index);                    
                    
                    upper_index+=1;                    
                    decimal_day_data = sat_data.data["decimal_day"][index+upper_index];
                    
                    
            for ind in new_indices:
                self.indices[key_filename]["decimal_day"].add(ind);

            # print;                
            # print key_filename;
            # print self.indices[key_filename];
            
            if (self.indices[key_filename]["decimal_day"]!=set()):            
                print min(sat_data.data["decimal_day"][ list(self.indices[key_filename]["decimal_day"])])
                print max(sat_data.data["decimal_day"][ list(self.indices[key_filename]["decimal_day"])])

    # This find the closest value to the event's value for the same keyname        
    def compare_with_zerocrossing(self,key_varnames):
        
        for sat_data in self.satellite_data.dataset[:]:
 
            # Create new set in the indices dictionary
            key_filename = sat_data.data_filename[:-12];

            for key_varname in key_varnames:
                
                self.key_varnames.add(key_varname);
                self.indices[key_filename][key_varname] = set();
               
                # Get the satellite's data
                sat_data_array = sat_data.data[key_varname];
                retain=False;

                # Get the sign of the first datapoint 
                del_data = (sat_data_array[0]-self.event.data[key_varname])  ;
                new_sign1 = np.sign(del_data);
                
                # Search for sign change in the difference
                for index,sat_data_element in enumerate(sat_data_array):

                    old_sign = new_sign1;
                    
                    del_data = (sat_data_element-self.event.data[key_varname])  ;
                    new_sign1 = np.sign(del_data);

                    # Check for crossing new year boundary in "decimal_day" data
                    if(key_varname=="decimal_day"):

                        startOfThisYear = datetime(year=self.event.date.year, month=1, day=1);
                        startOfNextYear = datetime(year=self.event.date.year+1, month=1, day=1);
                        days_in_year = float((startOfNextYear-startOfThisYear).days);       

                        del_data = (sat_data_element+days_in_year-self.event.data[key_varname]);
                        new_sign2 = np.sign(del_data);

                    else:
                        new_sign2 = new_sign1;
                        
                    if((new_sign1!=old_sign) and (new_sign2==new_sign1)):
                        retain=True;
                        
                    # Check for change of sign, ignore sign change at start of datafile
                    if(retain):
                        self.indices[key_filename][key_varname].add(index);
                        # var_sets[count].add(index);
                        retain=False;
               
                
                    
            # print sat_lat[indices];
            # print long_zero_crossings;
                    
        # self.print_indices();
        return 0;
    
    def print_indices(self):
        # Print indices         
        for key_filename in self.key_filenames:
            print;
            print key_filename;
            print self.indices[key_filename];
        
    def intersect_indices_set(self):
        # Save the set of indices of datapoints that are common to all search terms
        for key_filename in self.key_filenames:
            self.indices[key_filename]["intersection"] = self.indices[key_filename]["decimal_day"];
            for key_varname in self.key_varnames:
                self.indices[key_filename]["intersection"].intersection_update(self.indices[key_filename][key_varname]);
                print key_filename,key_varname, len(self.indices[key_filename]["intersection"])
        self.key_varnames.add("intersection");

    def compare_dates(self,tol):
        

        for sat_data in self.satellite_data.dataset[:]:
 
            # Create new set in the indices dictionary
            key_filename = sat_data.data_filename[:-12];
            self.indices[key_filename]['date'] = set();
            
            for variable_name in sat_data.data:
                self.selected_data[key_filename][variable_name] = list();

            # Search for sign change in the difference
            for index,sat_data_element in enumerate(sat_data.data['decimal_day']):

                # Find out how many days are in the year
                data_year = sat_data.data['year'][index];
                epoch = datetime(data_year - 1, 12, 31);
                result = epoch + timedelta(days=float(sat_data_element));

                # Date difference
                date_difference = result - self.event.date;    
                del_data = date_difference.seconds/(60.0*60.0*24.0);

                if(del_data<tol):
                    retain=True;
                        
                # Check for change of sign, ignore sign change at start of datafile
                if(retain):
                    retain=False;
                    self.indices[key_filename]['date'].add(index);
                    for variable_name in sat_data.data:
                        self.selected_data[key_filename][variable_name].append(sat_data.data[variable_name][index]);
                
        return 0;


   
    def compare_with_tolerance(self,key_varnames,tol):
        
        for sat_data in self.satellite_data.dataset[:]:
 
            # Create new set in the indices dictionary
            key_filename = sat_data.data_filename[:-12];

            for key_varname in key_varnames:
                
                self.key_varnames.add(key_varname);
                self.indices[key_filename][key_varname] = set();

                # Get the satellite's data
                sat_data_array = sat_data.data[key_varname];
                retain=False;
                
                # Search for sign change in the difference
                for index,sat_data_element in enumerate(sat_data_array):

                    del_data = math.fabs(sat_data_element-self.event.data[key_varname]);
                    
                    if(del_data<tol):
                        retain=True;
                        
                    # Check for change of sign, ignore sign change at start of datafile
                    if(retain):
                        self.indices[key_filename][key_varname].add(index);
                        retain=False;
                
        return 0;

    def compare(self,key_varnames,tol=0):
        if(tol==0):
            self.compare_with_zerocrossing(key_varnames);
        else:
            self.compare_with_tolerance(key_varnames,tol);



























