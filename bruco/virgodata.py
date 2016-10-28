# Brute force coherence (Gabriele Vajente, 2016-10-26)
# Virgo-centric functions (Bas Swinkels, Irene Fiori)

import sys
import numpy as np
# test to replace nap, hard-coded for now 
sys.path.append('/users/swinkels/deploy/PythonVirgoTools/trunk/src')
    
from virgotools import getChannel, FrameFile
  
# return the list of channels    
def get_channel_list(opt, minfs, gpsb):
    channels = []
    sample_rate = []
    channelMain = 'V1:' + opt.channel
    with FrameFile('raw') as ffl:
        with ffl.get_frame(gpsb) as frame:
            for adc in frame.iter_adc():
                fr = int(adc.contents.sampleRate)
                channel = str(adc.contents.name)
                channels.append(channel)	    
                sample_rate.append(fr)
      
    channels, sample_rate = (list(t) for t in zip(*sorted(zip(channels, sample_rate))))
    nchold = len(channels)  # number of channels before the exclusions
    #print "Found %d channels before the exclusions\n\n" % nchold
    return np.array(channels), np.array(sample_rate)

# Function to get data from raw files
def getRawData(ch, gps, dt):
    """Read data from RAW file:
    ch  = channel name
    gps = gps time
    dt  = duration
    """
    with getChannel('raw', ch, gps, dt) as x:
        return x.data, x.fsample
