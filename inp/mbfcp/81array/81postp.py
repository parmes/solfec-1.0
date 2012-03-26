
import sys
import csv
import sets
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl

# --- user parameters -----------------------------------------
argv = NON_SOLFEC_ARGV()

if argv == None:
  print '---------------------------------------------------------------------------------------'
  print 'SYNOPSIS: solfec 81postp.py path/to/file_1.thv label_1 path/to/file_2.thv label_2 [...]'
  print '---------------------------------------------------------------------------------------'
  print 'Paths and labels can be given in any combination, only their order matters.'
  print 'For example this is also fine: solfec 81postp.py path1 path2 label1 label2.'
  print 'You must first run analysis for 81array.py and then run it again print in'
  print 'read mode to extract the  *.thv file.'
  print '---------------------------------------------------------------------------------------'
  sys.exit ()

success = 0
vpath = []
vlabel = []
for str in argv:
  if str.endswith ('.thv'):
    try:
      f = open (str, 'r')
      success = 1
      f.close ()
      vpath.append (str)
    except: pass
  else: vlabel.append (str)

if not success:
  print 'ERROR: no *.thv file path has been passed'
  sys.exit ()

if len(vpath) != len(vlabel):
  print 'ERROR: number of input paths does not match number of labels'
  sys.exit ()

print 'Input files:', vpath
print 'Labels:', vlabel

#  --- per-analysis parameters --------------------------------
analysis = '81array' # analysis name
win_len = 0.1        # float - window length, in Hz, used for processing
swp_rate = 0.1       # float - input sweep rate, in Hz/s
start_freq = 3.0     # float - input start frequency, in Hz
# --------------------------------------------------------------

# --- configurable parameters, which shouldn't be changed ---
FBlabels = ['FB1(0)(0)',
'FB1(3)(4)', # B
'FB2(2)(3)', # E
'FB2(1)(2)', # C
'FB2(2)(2)', # T
'FB1(3)(3)', # F
'FB2(0)(2)', # G
'FB2(1)(1)', # J
'FB1(3)(1)', # Q
'FB2(1)(3)', # K
'FB1(3)(2)', # M
] # Fuel bricks results are required for

IBlabels = ['IB2(3)(2)', # 7
'IB2(1)(0)', # 18
'IB2(4)(2)', # 13
'IB2(4)(0)', # 17
] # Interstital bricks results are required for

boundarylabel = 'FB1(0)(0)' # label of a boundary brick - NOTE this MUST also be the first entry in 'FBlabels'
# --------------------------------------------------------------

__version__ = '1 DRAFT B1'

def WINDOWS(vallist, winlength):
  """ Returns a tuple-of-2tuples ( (start1, end1), (start2, end2), ..., (startX, endX))
      which defines the INDEXES of 'vallist' giving adjacent windows of length 'winlen'.
      Only windows which entirely fit in the vallist are returned.
  """
  
  idxs = []
  startidx = 0
  while startidx < len(vallist) - 1:
    endidx = startidx + winlength
    if endidx < len(vallist):
      idxs.append((startidx, endidx))
    startidx = endidx + 1  # note this will correctly fail the 'while' if endidx has been truncated
  return tuple(idxs)

def WINTIMES(times, windowindexes):
  """ Returns a list with the times of centre of windows in 'times' defined by 'windowindexes'.
      'windowindexes' should be produced using WINDOWS().
  """
  
  if not isinstance(windowindexes, tuple) and isinstance(windowindexes[1], tuple) and len(windowindexes[1]) == 2:
    raise TypeError, 'Second argument to WINTIMES must be a tuple-of-2tuples returned by WINDOWS()'
  
  ts = []
  for w in windowindexes:
    lowidx =  w[0]
    highidx = w[1]
    t = (times[lowidx] + times[highidx]) / 2.0
    ts.append(t)
  return ts
 
def PEAKWIN(series, windowindexes):
  """ Return a list where each entry is the maximum of absolute values in 'series' over the windows
      specified by the indexes in 'windowindexes'.
      'windowindexes' should be produced using WINDOWS()
  """
  
  if not isinstance(windowindexes, tuple) and isinstance(windowindexes[1], tuple) and len(windowindexes[1]) == 2:
    raise TypeError, 'Second argument to AVERWIN must be a tuple-of-2tuples returned by WINDOWS()'

  series = [abs(s) for s in series]
  vs = []
  for w in windowindexes:
    lowidx =  w[0]
    highidx = w[1]
    v = max(series[lowidx:highidx])
    vs.append(v)
  return vs

sys.stdout.flush()

# build results request
requests = []
for bl in FBlabels:

  requests.append( (bl, 'VX') )
  requests.append( (bl, 'VX') )

for bl in IBlabels:

  requests.append( (bl, 'VX') )
  requests.append( (bl, 'VX') )

# create figure handle
fig1 = plt.figure()

# configure plots
plotorder = ['FB2(2)(3)', # E
	     'FB2(2)(2)', # T
	     'IB2(1)(0)', # 18
	     'IB2(3)(2)', # 7
	     'FB2(1)(2)', # C
	     'FB1(3)(4)', # B
	    ] # Define plot order top-bottom, left-right

# read time series from file
for (path, lbl) in zip (vpath, vlabel):
  # output path
  outpath = path[0:len(path)-4]

  # read time series
  f = open (path, 'r')
  thv = pickle.load (f)
  f.close ()

  # times at which results are specified in database
  dbtimes = thv[0]    
  dwell = dbtimes[0] # length in seconds of constant-frequency dwell at start of analysis
  analysis_end = dbtimes [len(dbtimes)-1]

  # calculate window stuff:
  # (note no dwell correction is needed here as we've only read data from after the dwell)
  dbperiod = dbtimes[1] - dbtimes[0]                                   # period (s) between frames in the output database
  winsize = int((float(win_len) / float(swp_rate)) * (1.0 / dbperiod)) # number of data frames in each window
  windows = WINDOWS(dbtimes, winsize)                         # window start/end indices
  wintimes = WINTIMES(dbtimes, windows)                       # corresponding times of window centres
  winfreqs = [start_freq + t * swp_rate for t in wintimes]    # corresponding frequency at every window centre time

  # set up qa text
  qatext1 = 'analysis:%s, postprocess ver:%s' % (analysis, __version__)
  qatext2 = 'Output interval=%gs, window length=%gHz, entered sweep rate=%gHz/s => n=%ipoints' % (dbperiod, win_len, swp_rate, winsize)
  qatext3 = 'Data read from %gs to %gs' % (dwell, analysis_end)
  qatext = qatext1 + '\n' + qatext2 + '\n' + qatext3

  # open a file for text output
  with open(outpath + '_VX_UNVERIFIED.txt', 'w') as fout:
    fout.write(qatext + '\n')
    fout.write('\n')
    fout.write('\t'.join(['Body', 'Max(peak)', 'Freq']))
    fout.write('\n')
    
    # perform per-body calculations:
    bodyVs = {}      # 'raw' body centre velocity for each time in database
    peakVs = {}      # maximum-absolute values in windows
    nor_pkVs = {}    # normalised values series using peak in window
    
    for i in range(0, len(requests), 2):    # step through requests 2 at a time
      
      # get the body label for these 2 requests
      blabel = requests[i][0]
      
      # get time histories for both points
      tha = thv[i + 1] # first is +1 due to time being at results[0]
      thb = thv[i + 2]
      
      # calculate & store body centre velocities by averaging results for points
      bodyVs[blabel] = [(av + bv) / 2.0 for av, bv in zip(tha, thb)]
    
      # calculate series using functions over windows
      peakV = PEAKWIN(bodyVs[blabel], windows)
    
      # push these into dictionaries - need to do this now as we will access it immediately below
      peakVs[blabel] = peakV
    
      # normalise
      nor_pkV = [Vout / Vin for Vin, Vout in zip(peakVs[boundarylabel], peakVs[blabel])]
      
      # push these into dictionaries
      nor_pkVs[blabel] = nor_pkV
      
      # calculate max of series
      max_nor_pkV = max(nor_pkV)
      
      # calculate input frequency above max occurs at
      max_pk_idx = nor_pkV.index(max_nor_pkV)   # position in array of max
      max_pk_freq = winfreqs[max_pk_idx]          # input frequency of max
      
      # write max and input freq to file
      fout.write('%s\t%#g\t%g\n' % (blabel, max_nor_pkV, max_pk_freq))

    fout.write('End of data\n')
        
    ax = [] # will fill with axes objects for subplots as we generate them - NB is zero-based!
    for i, p in enumerate(plotorder):
      newax = fig1.add_subplot(3,2, i + 1) # row, col, plot_number (starts at 1, inconveniently, goes left-right top-bottom)
      newax.plot(winfreqs, nor_pkVs[p], '-', label=lbl)
      newax.set_title(p).set_fontsize('small')
      newax.grid(True)
      newax.set_ylim(0.0, 4.0)  # to match Fig 5 of C33/C34/PSD/213/220
      newax.set_xlim(0.0, 15.0) # ditto
      newax.legend(prop={'size':'small'})
      ax.append(newax)
    
    # add axes labels to edge subplots only
    ax[4].set_xlabel('Freq (Hz)', size='small')
    ax[5].set_xlabel('Freq (Hz)', size='small')
    ax[0].set_ylabel('Norm. velocity (-)', size='small')
    ax[2].set_ylabel('Norm. velocity (-)', size='small')
    ax[4].set_ylabel('Norm. velocity (-)', size='small')
    
    # add legend to first plot only
    
    # add titles and tidy output
    fig1.subplots_adjust(hspace=0.35, bottom=0.15, left=0.10, right=0.95)
    fig1.suptitle('%s: Normalised VX (boundary=%s) for %gHz windows' % (analysis, boundarylabel, win_len))
    fig1.text(0.01, 0.01, qatext).set_fontsize('small')
    
  # fout gets closed here by end of 'with' block

# plot experimental data 
for i, p in enumerate(plotorder):
  data = csv.reader (open ('inp/mbfcp/81array/' + p + '.csv', 'rb'))
  x = []
  y = []
  for row in data:
    x.append (row [0])
    y.append (row [1])
  newax = fig1.add_subplot(3,2, i + 1) # row, col, plot_number (starts at 1, inconveniently, goes left-right top-bottom)
  newax.plot(x, y, '-', label='EXP')
  newax.set_title(p).set_fontsize('small')
  newax.grid(True)
  newax.set_ylim(0.0, 4.0)  # to match Fig 5 of C33/C34/PSD/213/220
  newax.set_xlim(0.0, 15.0) # ditto
  newax.legend(prop={'size':'small'})
  ax.append(newax)

# find common output name
def makeset(path):
  val = []
  i = 1
  for j in path.split ('_'):
    if j.endswith ('thv'): continue
    val.append ((i, j))
    i = i + 1
  return sets.Set (val)

set1 = makeset (vpath [0])
for path in vpath:
  set2 = makeset (path)
  set1 = set1.intersection (set2)

list = []
for item in set1:
  list.append (item)
list.sort ()

plotpath = ''
for item in list:
  plotpath += item [1]

plotpath += '_VX_UNVERIFIED.eps'

# save figure
fig1.savefig(plotpath)

# record we've successfully saved this:
print 'Saved plot %s' % plotpath

# show figure
plt.show()
