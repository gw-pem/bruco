Brute force coherence (Gabriele Vajente, 2017-01-23 vajente@caltech.edu)

Command line arguments (with default values)

--ifo                                     interferometer prefix [H1, L1, V1] 
                                          (no default, must specify)

--channel=OAF-CAL_DARM_DQ                 name of the main channel

--excluded=bruco_excluded_channels.txt    file containing the list of channels excluded 
                                          from the coherence computation

--gpsb=1087975458                         starting time

--length=180                              amount of data to use (in seconds)

--outfs=8192                              sampling frequency of the output results 
                                          (coherence will be computed up to outfs/2 
                                          if possible)

--minfs=512                               skip all channels with samplig frequency 
                                          smaller than this

--naver=100                               number of averages to compute the coherence

--dir=bruco                               output directory

--top=100                                 for each frequency, save to cohtab.txt and 
                                          idxtab.txt this maximum number of coherence 
                                          channels

--webtop=20                               show this number of coherence channels per 
                                          frequency, in the web page summary

--plot=png                                plot format (png, pdf, none)

--nproc                                   number of processes to use (if not specified,
                                          use all CPUs)

--calib                                   name of a text file containing the calibration 
                                          transfer function to be applied to the target 
                                          channel spectrum, in a two column format 
                                          (frequency, absolute value)

--xlim                                    limits for the frequency axis in plots, in the 
                                          format fmin:fmax

--ylim                                    limits for the y axis in PSD plots, in the 
                                          format ymin:ymax

--tmp=~/tmp                               temporary file directory where cache files 
                                          will be saved if needed

Example:
./bruco.py --ifo=H1 --channel=CAL-DELTAL_EXTERNAL_DQ 
           --calib=share/lho_cal_deltal_calibration.txt 
           --gpsb=1111224016 --length=600 --outfs=4096 --naver=100 
           --dir=./bruco_1111224016 --top=100 --webtop=20 --xlim=1:1000 
           --ylim=1e-20:1e-14 --excluded=share/bruco_excluded_channels.txt"""

CHANGELOG:

2015-01-29 - added linear detrending in PSD and coherence to improve low frequency bins
2015-02-16 - split coherence computation into PSDs and CSDs to remove redundant 
             computation of main channel PSD (still room to improve here)
           - removed upsampling of channels with lower sampling rate (now coherences are
             always computed with the minimum possible number of samples)
2015-02-18 - further splitting of coherence into primitive FFTs
2015-02-19 - allow selection of plot output format (PDF or PNG)
2015-02-24 - corrected typo in options lenght -> length
           - implemented parallel processing
2015-06-11 - using gw_data_find to locate the GWF files
2015-06-16 - added calibration transfer function option
           - added expanduser to all paths to allow use of ~
2015-08-19 - now bruco removes the temporary files that contains the list of channels
           - timing information file is saved into the output directory
2015-01-01 - x3 faster plotting based on Michele Valentini's code
2016-09-19 - fixed bug in output table
2016-10-18 - added command line options for excluded channels file and temporary folder, 
             updated default values   
2016-10-26 - added portability to Virgo environment
           - removed old timing code
           - when invoked without arguments, print help
2016-10-28 - reverted to non packaged structure for simplicity
2016-11-11 - using direct frame file access to read data, instead of LAL cache (a couple 
             of minutes faster)
2017-01-04 - using resample function if sampling frequency ratio is not integer
2017-01-05 - explicitly removing main channel from aux channel list
2017-01-23 - if the decimation ratio is greater than 10, decimate in mutiple steps to
	     avoid numerical instabilities
2017-03-08 - Virgo sampling rate is not always a power of 2, and this caused a crash
             for some auxiliary channels with sampling rate lower than the desired
             output. Implemented an upsamping of the aux channel to solve the issue
2017-04-04 - mask frequencies above Nyquist

