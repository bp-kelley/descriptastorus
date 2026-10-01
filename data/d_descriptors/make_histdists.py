"""Make histograms from all the distributions in the data directory
n.b. only make new distributions

Run from this directory:  python make_histdists.py
New histograms are written into descriptastorus/descriptors/hists.py
"""
import os
import sys
import gzip
import numpy
from numpy import inf, nan

# Load hists.py by path so this works without installing descriptastorus
# (importing the package pulls in pandas, pandas_flavor, rdkit ...)
HISTS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "descriptastorus", "descriptors",
                          "hists.py")
HISTS_FILE = os.path.normpath(HISTS_FILE)
_ns = {}
with open(HISTS_FILE) as f:
    exec(f.read(), _ns)

changed = False
histdists = _ns['hists']
for fname in sorted(os.listdir('.')):
    head, ext = os.path.splitext(fname)
    if ext == ".gz":
        name = head.replace("d_", "", 1)
        if name in histdists:
            print(f"Skipping {name}", file=sys.stderr)
            continue

        with gzip.open(fname) as f:
            txt = f.read()
        dist = numpy.asarray(eval(txt), dtype=float)
        # drop inf/-inf/nan, numpy.histogram can't bin them
        dist = dist[numpy.isfinite(dist)]
        n = min(1000, len(set(dist)))
        hist, xaxis = numpy.histogram(dist, bins=n)
        total = hist.sum()
        bins = []
        last = 0.0
        for value, x in zip(hist, xaxis):
            assert value >= 0
            last += value
            bins.append((float(x), float(last / total)))

        print(f"Adding {name}", file=sys.stderr)
        histdists[name] = bins
        changed = True

if changed:
    print(f"Writing updated histograms to {HISTS_FILE}", file=sys.stderr)
    text = repr(histdists)
    text = text.replace("],", "],\n\t")
    with open(HISTS_FILE, 'w') as f:
        f.write(f"hists = {text}")
else:
    print("No new datafiles added", file=sys.stderr)
