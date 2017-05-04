.. |Logo| image:: https://bitbucket.org/charade/zoomx/raw/master/doc/images/logo.png
   :alt: logo.png
   :height: 50px
   :width: 50px

.. |Rationale| image:: https://bitbucket.org/charade/zoomx/raw/master/doc/images/rationale.png
   :alt: rationale.png
   :height: 450px
   :width: 540px

.. |Pipeline| image:: https://bitbucket.org/charade/zoomx/raw/master/doc/images/pipeline.png
   :alt: pipeline.png
   :height: 320px
   :width: 520px

|Logo| ZoomX - Single Molecule Based Rearrangement Analysis with Linked Read Sequencing
==========================================================================================

.. contents:: Navigator

QUICK LINKS
-----------

`Examples <Example>`__

`Manuals <Manual>`__

`FAQ <FAQ>`__

INTRODUCTION
------------
    
  ZoomX (Single Molecule Based Rearrangement Analysis with Linked Read Sequencing)

  |Rationale|

  Figure 1. The principle and workflow of ZoomX.
  
  |Pipeline|

  Figure 2. The bioinformatics pipeline of ZoomX.
    
  Currently the package works for Linux (tested with Ubuntu). 
  It might also work for Mac with Brew/Ports or Windows with Cygwin (not tested).
  ZoomX documentation is on wiki and it is available:
  http://bitbucket.org/charade/zoomx/wiki

TRYOUT
------------

  *AWS*:

    ZoomX is available for tryout in a live Amazon EC2 instance.
    Simply log on with SSH and find the package located in ~/setup/zoomx .
    ssh bkd@ec2-54-190-36-93.us-west-2.compute.amazonaws.com # secret is the same as user
    Testing data and scripts are available in ~/setup/zoomx/test .

  *VirtualMachine*:

    The Amazon Machine Image (ami) of the same EC2 instance is downloadable from Amazon S3 bucket.
    https://s3-us-west-2.amazonaws.com/lixiabucket/zoomxDemo.ova

INSTALL
------------
    
  *Dependencies*:

    Python2(>=2.7)  http://www.python.org/

      pybedtools, pysam, pandas, tables, scipy

    R(>=3.3)  http://r-project.org
      
      optparse, GenmicRanges, ggbio, rjson, RColorBrewer 
      
    Samtools (>0.19) http://samtools.org
      executables have to be in $PATH and within the sve_py virtual env
      there is samtools api comes with Pysam, however it has to be able to compile on your system with dafault setting. otherwise expect the "csamtools.c" not found error. e.g. libcurse is required for build samtools but often missing in some Linux.
      
    Bedtools (>2.20) http://bedtools.org

    An example install dependency install script for Ubuntu 16 is here:

      https://bitbucket.org/charade/zoomx/raw/master/test/ubuntu.install
    
  *InstallationFree*:

    ZoomX can be used without installation as long as required commands are all in $PATH.
    see manual for more informations

  *Installation*
    
    virtualenv based non-root installation with python site separation. 
    virtualenv command is standard with Python 2.7 or later. 
    If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. Possibly as simple as:
    
    .. code::

      sudo easy_install pip
      sudo pip install virtualenv

    Ask your IT manager to help install it for you if you have permission difficulties.
    When your system python has virtualenv, make sure your $PYTHONPATH is set to empty. 
    Follow steps below:
    
    .. code::
      
      sudo pip install --upgrade pip
      virtualenv zoomx_vpy

    Now you can activate this virtual python for installing ZoomX python dependencies:
    
    .. code::
      
      # in zoomx_vpy>
      source zoomx_vpy/bin/activate 
      pip install scipy
      pip install tables
      pip install pandas
      pip install pysam
      pip install pybedtools

    Optionally, on most platforms, the python dependencies will be automatically setup by ZoomX:
    
    .. code::

      # in zoomx_vpy>
      git clone git@bitbucket.org:charade/zoomx.git
      cd zoomx
      git submodule update --init --recursive 
      python setup.py install 

    Now the ZoomX executables will be available from "$PWD/zoomx_vpy/bin". 
    Because you installed ZoomX via virtualenv, remember to activate the virtualenv every time you use ZoomX. 
    Also export the environmental variable $ZOOMX_BIN=$PWD/zoomx_vpy/bin and add it to your $PATH

EXECUTABLES      
-----------

  parse10x      --  parsing 10X genomics (Chromium) bam
    
  gridScan      --  scalable Poisson grid whole genome scan 
    
  refineBreak   --  refine and report junction candidates 

  filterMolecule.R -- filter molecules based on coverage (requires R>3.3 and optparse, IRanges, GenomicRanges from CRAN and BioConductor)

  plotMolecule.R  --  plot molecules supporting a junction (requires R>3.3 and optparse, rjson, ggplot2, RColorBrewer, GenomicRanges, ggbio from CRAN and BioConductor)

  getMolecule.sh  --  output molecules supporting a junction 
 
  bedpe2loupe.sh  --  convert bedpe coordinates to loupe coordinates for visualization 

USAGE
------------
    
  By default all above executables will be available from $ZOOMX_BIN/ .
  Use '-h' to read script-wise usage.
  See Wiki pages for more details.

WIKI
-------------

  ZoomX's https://bitbucket.org/charade/zoomx/wiki page is a growing resource for manuals, FAQs and other information. This is a MUST read place before you actually using the ZoomX tool. These documentations are also openly editable in reSructuredText format. You are more than welcome to contribute to this ongoing documentation.
    
FAQ
-------------

  http://bitbucket.org/charade/zoomx/wiki/FAQ
    
BUG
-------------

  https://bitbucket.org/charade/zoomx/issues

CONTACTS
-------------

  Questions and comments shall be addressed to lixia at stanford dot edu.

CITATIONS
---------

  The ZoomX manuscript is under review.
