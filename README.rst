README
--------

INTRODUCTION
============
    
    ZoomX(Rearrangement Analysis with 10X Genomics Linked Read Sequencing)
    
    Currently the package works for Linux (tested with Ubuntu). 
    It might also work for Mac with Brew/Ports or Windows with Cygwin (not tested).
    ZoomX documentation is on wiki and it is available:
    http://bitbucket.org/charade/zoomx/wiki

INSTALL
========
    
    *Dependencies* (and that of):
    
    ::

      Python2(>=2.7)  http://www.python.org/

      Pysam  http://pysam.readthedocs.org/en/latest/
      
        Pysam prerequisites:
        Samtools (>0.19) http://samtools.org
          (executables have to be in $PATH and within the sve_py virtual env)
          (there is samtools api comes with Pysam, however it has to be able to compile on your system with dafault setting. otherwise expect the "csamtools.c" not found error. e.g. libcurse is required for build samtools but often missing in some Linux.)
      
      Bedtools http://bedtools.org 
      

    *Installation*
    
    virtualenv based non-root installation and python site separation. 
    
    (1) Install Engine through system/site python and virtualenv
      
      (1.1) virtualenv command is standard with Python 2.7 or later. If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. Possibly as simple as:
      
      ::

        sudo easy_install pip
        sudo pip install virtualenv

      Ask your IT manager to help install it for you if you have permission difficulties.
      
      (1.2) When your system python has virtualenv, make sure your $PYTHONPATH is set to empty and follow steps below:
      
      ::
        
        >virtualenv zoomx_vpy --no-site-packages

      (1.3) Then you can activate this virtual python for installing SVEngine:
      
      ::
        
        >source zoomx_vpy/bin/activate

      (1.4) Now under your virtualenv, the dependencies will be automatically setup by ZoomX installation itself:
      
      ::
        
        zoomx_vpy> python setup.py install 

      (1.5) On some platforms, you might run into error and it prompts you to install Cython, just install Cython and try again:
      
      ::
        
        zoomx_vpy>pip install -U cython
      
      (1.6) Now the SVEngine executables will be available from "$PWD/zoomx_vpy/bin". Because you installed SVEngine via virtualenv, remember to activate the virtualenv first every time you use SVEngine. Also export the environmental variable $ZOOMX_BIN=$PWD/zoomx_vpy/bin
       

EXECUTABLES      
===========

    parse10x      --  parsing 10X genomics (Chromium) bam
    
    gridScan      --  scalable Poisson grid whole genome scan 
    
    refineBreak   --  refine and report junction candidates 

    getMolecule   --  output molecules supporting a junction 

    plotMolecule  --  plot molecules supporting a junction 


USAGE
========
    
    (1) By default all above executables will be available from $ZOOMX_bin/ .
    Use '-h' to read script-wise usage.

    (2) Do a Sanity check for installation

    ::

      cd test
      ./test_mf.sh
    
WIKI
========
    http://bitbucket.org/charade/zoomx/wiki/Home
    
FAQ
========
    http://bitbucket.org/charade/zoomx/wiki/FAQ
    
BUG
========
    https://bitbucket.org/charade/zoomx/issues

CONTACT
========
    lixia at stanford dot edu