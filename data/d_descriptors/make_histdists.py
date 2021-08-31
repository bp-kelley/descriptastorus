"""Make histograms from all the distributions in the data directory
n.b. only make new distributions
"""
import os
import gzip
import numpy
import functools, operator
from numpy import inf,nan
import bisect
from descriptastorus.descriptors import hists

def histcdf(x, dist, N):
    p = bisect.bisect(dist, (x,))
    if p < N:
        return dist[p][1]
    return 1.0

changed = False
histdists = hists.hists
for f in os.listdir('.'):
    head, ext = os.path.splitext(f)
    if ext == ".gz":
        name = head.replace("d_", "")
        if name in histdists:
            print(f"Skipping {name}", file=sys.stderr)
            continue
        
        with gzip.open(f) as f:
            txt = f.read()
            dist = eval(txt)
            if inf in dist:
                dist = [x for x in dist if numpy.isfinite(x)]
            n = min(1000,len(set(dist)))
            hist, xaxis = numpy.histogram(dist, bins=n)
            total = functools.reduce(operator.add, hist)
            bins = []
            last = 0.0
            for value, x in zip(hist, xaxis):
                assert value >= 0
                last += value
                bins.append((x, last/total))
            N = len(dist)
            for v in dist:
                histcdf(v, dist, N)
                
            histdist[name] = bins
            changed = True

if changed:            
    print("Writing temporary histsdist to data directory", file=sys.stderr)
    text = repr(histdist)
    text = text.replace("],", "],\n\t")
    filename = hists.__file__
    if os.path.splitext(filename) == ".pyc":
        filename = filename[:-1]
        
    open(filename, 'w').write(f"hists = {text}")
else:
    print("No new datafiles added", file=sys.stderr)
