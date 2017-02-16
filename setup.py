#!/usr/bin/env python

"""ZoomX Package.
"""

from setuptools import setup, find_packages
import os, sys, subprocess

doclines=__doc__.splitlines()

generateVersion="""echo 'def run():\n\tprint \"'$(cat VERSION.txt)$(git log --pretty=format:'%h' | head -n 1 | sed s/^/\(/g | sed s/$/\)/g)'\" \n\nif __name__ == \"__main__\":\n\trun()' >zoomx/zoomxVersion.py"""
subprocess.call(generateVersion, shell=True)

"""samtools"""

print "cd samtools && make"
os.system("cd samtools && make")
print "cp samtools/samtools %s" % sys.prefix
os.system("cp samtools/samtools %s/bin" % sys.prefix)

"""bedtools"""

print "cd bedtools2 && make"
os.system("cd bedtools2 && make")
print "cp bedtools2/bin/* %s/bin" % sys.prefix
os.system("cp bedtools2/bin/* %s/bin" % sys.prefix)

#"""scripts"""
#print "cp scripts/* %s/bin" % sys.prefix
#os.system("cp scripts/* %s/bin" % sys.prefix)

dist=setup(name="zoomx",
    version="1.0.0",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Li Xia",
    author_email="lixia@stanford.edu",
    url="https://bitbuket.com/zoomx",
    license="TBD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    zip_safe=False,
    install_requires=["python >= 2.7","pandas >= 0.16", "pysam >= 0.7",],
    provides=['zoomx'],
    py_modules = ['zoomx.parse10x', 'zoomx.gridScan', 'zoomx.refindBreak'],
		scripts = ['scripts/getMolecule.sh', 'scripts/plotMolecule.R', \
							 'scripts/bedpe2loupe.sh', 'scripts/getSomatic.sh', \
               'scripts/plot10x.R', 'scripts/plotCircus.R', \
               'scripts/plotCumCov.R', 'scripts/plotDenCov.R'],
    data_files = [('',['README.rst','LICENSE','VERSION.txt'])],
    entry_points = { 
      'console_scripts': [
				'parse10x = zoomx.parse10x:run',
				'gridScan = zoomx.gridScan:run',
				'refindBreak = zoomx.refineBreak:run',
				'zoomxVersion = zoomx.zoomxVersion:run',
      ]
   },
)
