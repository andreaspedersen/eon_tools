import os, sys

thisdir, f = os.path.split(os.path.realpath(__file__))
libdir = os.path.join(thisdir, 'lib')

if libdir not in sys.path:
    sys.path.insert(0, libdir)
