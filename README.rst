.. |Logo| image:: https://bitbucket.org/charade/zoomx/src/master/doc/images/logo.png?fileviewer=file-view-default
   :alt: logo.png
   :height: 50px
   :width: 50px

.. |Rationale| image:: https://bitbucket.org/charade/zoomx/src/master/doc/images/rationale.png?fileviewer=file-view-default
   :alt: rationale.png
   :height: 450px
   :width: 540px

.. |Pipeline| image:: https://bitbucket.org/charade/zoomx/src/master/doc/images/pipeline.png?fileviewer=file-view-default
   :alt: pipeline.png
   :height: 320px
   :width: 520px

|Logo| ZoomX - Single Molecule Based Rearrangement Analysis with Linked Read Sequencing
==========================================================================================

.. contents:: Navigator

QUICK LINKS
-----------

`Examples <Example>`_

`Manuals <Manual>`_

`FAQ <FAQ>`_

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
    
  *Dependencies* (and that of):

    Python2(>=2.7)  http://www.python.org/

    Pysam  http://pysam.readthedocs.org/en/latest/
      
    Samtools (>0.19) http://samtools.org
      (executables have to be in $PATH and within the sve_py virtual env)
      (there is samtools api comes with Pysam, however it has to be able to compile on your system with dafault setting. otherwise expect the "csamtools.c" not found error. e.g. libcurse is required for build samtools but often missing in some Linux.)
      
    Bedtools http://bedtools.org       
    
    R  http://r-project.org

  *NoInstallation*:

    ZoomX can be used without installation as long as required commands are all in $PATH.
    see manual for more informations

  *Installation*
    
    virtualenv based non-root installation and python site separation. 
    
    (1) Install through system/site python and virtualenv
      
      (1.1) virtualenv command is standard with Python 2.7 or later. If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. Possibly as simple as:
      
      .. code::

        sudo easy_install pip
        sudo pip install virtualenv

      Ask your IT manager to help install it for you if you have permission difficulties.
      
      (1.2) When your system python has virtualenv, make sure your $PYTHONPATH is set to empty and follow steps below:
      
      .. code::
        
        sudo pip install --upgrade pip
        virtualenv zoomx_vpy

      (1.3) Then you can activate this virtual python for installing ZoomX:
      
      .. code::
        
        # in zoomx_py>
        source zoomx_vpy/bin/activate 
        pip install scipy
        pip install tables
        pip install pandas
        pip install pysam
        pip install pybedtools

      (1.4) Optionally, under your virtualenv, the dependencies will be automatically setup by ZoomX installation itself:
      
      .. code::

        # in zoomx_py>
        git clone git@bitbucket.org:charade/zoomx.git
        cd zoomx
        git submodule update --init --recursive 
        python setup.py install 

      (1.5) Now the ZoomX executables will be available from "$PWD/zoomx_vpy/bin". Because you installed ZoomX via virtualenv, remember to activate the virtualenv first every time you use ZoomX. Also export the environmental variable $ZOOMX_BIN=$PWD/zoomx_vpy/bin
       

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
    
  (1) By default all above executables will be available from $ZOOMX_bin/ .
    Use '-h' to read script-wise usage.

  (2) Do a Sanity check for installation

  .. code::

     cd test
     ./simple Your/Path/To/ZoomX
    
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
