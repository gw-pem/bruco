# Brute force coherence (Gabriele Vajente, 2016-10-18)
# LIGO-centric functions

#from glue import lal
import os
from pylal import frutils, Fr
import subprocess
from pylab import *
import numpy

# global variables for LIGO data I/O
files = []

# wrapper around the LIGO function to find where data is, returns a list of files
# def find_data(observatory, gpsb, gpse):
#     o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
#                           "-t", observatory[0] + "1_R", "-s", str(gpsb),
#                           "-e", str(gpse), "-u", "file"],
#                           stdout=subprocess.PIPE).communicate()[0]
#     return o.splitlines()

# wrapper around the LIGO function to find where data is, returns a list of files
def find_data_path(observatory, gpsb, gpse):
    pref = '/frame0/full'
    flen = int(32)
    gstt = int(int(gpsb / flen) * flen)
    gstp = int(int((gpse-0.5) / flen) * flen)
    fpaths = []
    for gpst in range(int(gstt), int(gstp+flen), int(flen)):
        print(str(gpst) + '\n')
        gdir = int(round(gpst / 100000))
        fpaths.append(pref + '/' + str(gdir) + '/K-K1_C-' + str(gpst) + '-' + str(flen) + '.gwf')
        
    return fpaths
    # o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
    #                       "-t", observatory[0] + "1_R", "-s", str(gpsb),
    #                       "-e", str(gpse), "-u", "file"],
    #                       stdout=subprocess.PIPE).communicate()[0]
    # return map(lambda x: x[16:], o.splitlines())

# get the list of channels
def get_channel_list(opt, gpsb):
    
    files0 = find_data_path(opt.ifo, gpsb, gpsb+1) 
    os.system('/users/tyamamoto/apps/deb-8.5/libframe-8.21/bin/FrChannels ' + files0[0] + ' > bruco.channels')
    f = open('bruco.channels')
    lines = f.readlines()
    channels = []
    sample_rate = []
    for l in lines:
        ll = l.split()
        if ll[0][1] != '0':
            # remove all L0/H0 channels
            channels.append(ll[0])
            sample_rate.append(int(ll[1]))
    channels = array(channels)
    sample_rate = array(sample_rate)
    # remove temporary channel list
    os.system('rm bruco.channels')
    
    return channels, sample_rate
    
# Function to get data from raw files
#def getRawData(channel, gps, dt):
#    """Read data from RAW file:
#    ch  = channel name
#    gps = gps time
#    dt  = duration
#    """
#    global ligo_c
#    global ligo_d
#    
#    buffer = ligo_d.fetch(channel, gps, gps+dt)
#    ch1 = numpy.array(buffer)
#    fs1 = len(ch1) / dt
#    
#    return ch1, fs1

def getRawData(channel, gps, dt):
    """Read data from RAW file:
    ch  = channel name
    gps = gps time
    dt  = duration
    """
	
    global files

    # get the list of frame files
    obs = channel[0:2]
    if len(files) == 0:
    	files = find_data_path(obs, gps, gps+dt)
    # read data from all files
    data = array([]) 
    for f in files:
	# get the frame file start GPS and duration from the name 
    	gps0 = int(f.split('-')[-2])
	gps1 = gps0 + int(f.split('-')[-1].split('.')[-2])
	# find the right time segment to load from this file
        gps0 = max(gps0, gps)
	gps1 = min(gps1, gps + dt)
	# read data and append
	x = Fr.frgetvect(f, channel, gps0, gps1-gps0)
	data = concatenate([data, x[0]])
    return data, int(1/x[3][0])

