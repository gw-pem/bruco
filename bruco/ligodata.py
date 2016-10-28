# Brute force coherence (Gabriele Vajente, 2016-10-18)
# LIGO-centric functions

from glue import lal
import os
from pylal import frutils, Fr
import subprocess
from pylab import *
import numpy

# global variables for LIGO data I/O
ligo_c = 0
ligo_d = 0

# wrapper around the LIGO function to find where data is, returns a list of files
def find_data(observatory, gpsb, gpse):
    o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
                          "-t", observatory[0] + "1_R", "-s", str(gpsb),
                          "-e", str(gpse), "-u", "file"],
                          stdout=subprocess.PIPE).communicate()[0]
    return o.splitlines()

# get the list of channels
def get_channel_list(opt, minfs, gpsb):
    global ligo_c
    global ligo_d
    
    # create data access objects if needed
    if ligo_c == 0:
        print ">>>>> Creating cache...."
        # Create the scratch directory if not exisitng
        try:
            os.stat(os.path.expanduser(opt.scratchdir))
            new_scratch = False
        except:
            os.mkdir(os.path.expanduser(opt.scratchdir))
            new_scratch = True
        # find the location of the data files
        files = find_data(opt.ifo, gpsb, gpsb+int(opt.dt))
        # create LAL cache
        ligo_c = lal.Cache.from_urls(files)
        ligo_d = frutils.FrameCache(ligo_c, scratchdir=os.path.expanduser(opt.scratchdir), 
                                        verbose=True)

    # read the list of channels and sampling rates from the first file
    firstfile = ligo_c[0].path
    os.system('/usr/bin/FrChannels ' + firstfile + ' > bruco.channels')
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
def getRawData(channel, gps, dt):
    """Read data from RAW file:
    ch  = channel name
    gps = gps time
    dt  = duration
    """
    global ligo_c
    global ligo_d
    
    # create data access objects if needed
    if ligo_c == 0:
        print ">>>>> Creating cache...."
        # Create the scratch directory if not exisitng
        try:
            os.stat(os.path.expanduser(opt.scratchdir))
            new_scratch = False
        except:
            os.mkdir(os.path.expanduser(opt.scratchdir))
            new_scratch = True
        # find the location of the data files
        files = find_data(opt.ifo, gpsb, gpsb+dt)
        # create LAL cache
        ligo_c = lal.Cache.from_urls(files)
        ligo_d = frutils.FrameCache(ligo_c, scratchdir=os.path.expanduser(opt.scratchdir), 
                                        verbose=True)

    buffer = ligo_d.fetch(channel, gps, gps+dt)
    ch1 = numpy.array(buffer)
    fs1 = len(ch1) / dt
    
    return ch1, fs1