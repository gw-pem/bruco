# Brute force coherence (Gabriele Vajente, 2015-06-09)

import nds2
import numpy
import os
from pylab import *
import time
from scipy.signal import cheby1, lfilter
import subprocess

# Compute the windowed FFTs of the channel x, each with npoints and overlapped by noverlap.
# Returns all the FFTs, not only the average.
def computeFFTs(x, npoints, noverlap, fs):
    """Compute the windowed FFTs of the channel x, each with npoints and overlapped by noverlap."""

    # compute the starting indexes of all segments
    step = npoints - noverlap
    ind = arange(0, len(x) - npoints + 1, step)
    nsegs = len(ind)

    # pre-compute window
    wind = numpy.hanning(npoints)

    # pre-allocate output data
    segs = zeros((npoints/2+1, nsegs), dtype=numpy.cfloat)

    # compute all FFTs
    for i in range(nsegs):
        segs[:,i] = rfft( detrend_linear(x[ind[i]:ind[i]+npoints]) * wind)

    # normalization to get PSD with |fft|^2
    segs = segs * sqrt(2.0/fs / sum(wind**2))

    # that's all
    return segs

# Convert GPS time to string time (not sure if the conversion constant is still good)
def gps2str(gps, zone='utc'):
    """Convert GPS to string time"""
    if (zone == 'utc'):
        tz = time.mktime(time.gmtime()) - time.mktime(time.localtime())
        return time.ctime(gps + 315964785 + tz) + ' UTC'
    else:
        return time.ctime(gps + 315964785)

# return the actual GPS time  (not sure if the conversion constant is still good)
def gps_now():
    """Return the actual GPS time."""
    return int(time.time() - 315964785)

# Build a color code from white to red depending on the value of the argument
def cohe_color(c):
    """Return an hex color code depending continuosly
    on input: 1 = red, 0 = white
    """
    if c == 0:
        return '#ffffff'
    else:
        if c == 1:
            return '#ff0000'
        else:
            s = hex(int((1.0 - c) * 256)).split('x')[1]
            if len(s) == 1:
                s = '0' + s
            return '#ff' + s + s

# Custom decimation function, copied a long time ago from somewhere on the web
def decimate(x, q, n=None, ftype='iir', axis=-1):
    """downsample the signal x by an integer factor q, using an order n filter

    By default, an order 8 Chebyshev type I filter is used or a 30 point FIR
    filter with hamming window if ftype is 'fir'.

    (port to python of the GNU Octave function decimate.)

    Inputs:
        x -- the signal to be downsampled (N-dimensional array)
        q -- the downsampling factor
        n -- order of the filter (1 less than the length of the filter for a
             'fir' filter)
        ftype -- type of the filter; can be 'iir' or 'fir'
        axis -- the axis along which the filter should be applied

    Outputs:
        y -- the downsampled signal

    """

    if type(q) != type(1):
        raise Error, "q should be an integer"

    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 4
    if ftype == 'fir':
        b = firwin(n+1, 1./q, window='hamming')
        y = lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = cheby1(n, 0.05, 0.8/q)

        y = lfilter(b, a, x, axis=axis)

    return y.swapaxes(0,axis)[::q].swapaxes(0,axis)

# This is a function to split a channel name in multiple lines if too long
def newline_name(s):
    if len(s) > 10:
        N = len(s)/10
        idx = []
        for i in range(N):
            try:
                idx.append(s.index('_', 10*(i+1)))
            except:
                pass
        if len(idx) !=0:
            newstr = ''
            for i in range(len(idx)):
                if i == 0:
                    newstr = s[0:idx[0]]
                else:
                    newstr = newstr + "<br>" + s[idx[i-1]:idx[i]]
            newstr = newstr + '<br>' + s[idx[-1]:]
            return newstr
        else:
            return s
    else:
        return s

# wrapper around the LIGO function to find where data is, returns a list of files
def find_LIGO_data(observatory, gpsb, gpse):
    o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
                          "-t", observatory[0] + "1_R", "-s", str(gpsb),
                          "-e", str(gpse), "-u", "file"],
                          stdout=subprocess.PIPE).communicate()[0]
    return o.splitlines()
