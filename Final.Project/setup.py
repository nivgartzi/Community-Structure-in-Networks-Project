
from setuptools import setup, Extension

setup(name='mykmeanssp',
      version='1.0',
      description='mykmeanssp setup.py',
      ext_modules=[Extension('mykmeanssp', sources=['spkmeans.c', 'spkmeansmodule.c'])])

# the extension name is like the name in the signaure: PyInit_XXX(void)   (a function in spkmeansmoudle.c)
