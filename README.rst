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

`Examples <https://bitbucket.org/charade/zoomx/wiki/Example>`__

`Manuals <https://bitbucket.org/charade/zoomx/wiki/Manual>`__

`FAQ <https://bitbucket.org/charade/zoomx/wiki/FAQ>`__

INTRODUCTION
------------
    
  ZoomX (Single Molecule Based Rearrangement Analysis with Linked Read Sequencing)

  - ZoomX analyzes single molecule whole genome sequencing data sets to identify long-range rearrangements.

  - ZoomX currently supports linked read sequencing data from 10X Genomics Chromium preparation.

  - ZoomX takes as input the BAM file output of LongRanger aligner and outputs a list of genomic junction sites in BEDPE format.

  - ZoomX's modeling and pipeline are illustrated in Figures 1 and 2.

  Figure 1. The principle and workflow of ZoomX.

  |Rationale|

  Figure 2. The bioinformatics pipeline of ZoomX.
  
  |Pipeline|
    
  Currently the package works for Linux (tested with Ubuntu). 
  It might also work for Mac with Brew/Ports or Windows with Cygwin (not tested).
  ZoomX documentation is on wiki and it is available:
  http://bitbucket.org/charade/zoomx/wiki

TRYOUT
------------

  An AWS EC2 instance and a preconfigured virtual machine are available for trying out ZoomX without installation efforts.

  *AWS*:

    ZoomX is available for tryout in a live Amazon EC2 instance.
    Simply log on with SSH and find the package located in ~/setup/zoomx .
    
    .. code:
    
      ssh bkd@ec2-54-190-36-93.us-west-2.compute.amazonaws.com # secret is the same as user
    
    Testing data and scripts are available in ~/setup/zoomx/test . Run ./clean to clean up the directory first. Examples are offered as germline3.run and somatic7.run, as explained in `Examples <https://bitbucket.org/charade/zoomx/wiki/Example>`__ .

    .. code:
   
      cd  ~/setup/zoomx/test 
      ./clean

  *VirtualMachine*:

    The Amazon Machine Image (ami) of the above EC2 instance is downloadable from Amazon S3 bucket.
    https://s3-us-west-2.amazonaws.com/lixiabucket/zoomxDemo.ova

INSTALL
------------
    
  *Dependencies*:

    Python2(>=2.7)  http://www.python.org/

      pybedtools, pysam, pandas, tables, scipy
      However, to install pysam, that api has also to be able to compile on your system with dafault setting. 
      Otherwise expect the "csamtools.c" not found error, sometimes because of libcurse is required for build samtools but often missing in Linux distros.

    R(>=3.3)  http://r-project.org
      
      CRAN: optparse, rjson, RColorBrewer
 
      Bioconductor: GenomicRanges, ggbio, GenomeInfoDb

      NOTE: using non-essential R scripts not listed in the EXECUTABLES section of this README file may incurr more dependencies, such as rhdf5, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg19, cowplot. 
      
    Samtools (>0.19) http://samtools.org
      samtools executables have to be in $PATH and within the sve_py virtual env
      NOTE: this is independent of the samtools api comes with Pysam. 
      
      
    Bedtools (>2.20) http://bedtools.org
      bedtools executables have to be in $PATH and within the sve_py virtual env
      
    An example dependency install script for Ubuntu is here:

      https://bitbucket.org/charade/zoomx/raw/master/test/ubuntu.install   

    
  *InstallationFree*:

    ZoomX can be used without installation as long as required commands are all in $PATH.
    see `Examples <https://bitbucket.org/charade/zoomx/wiki/Example>`__ and `Manuals <https://bitbucket.org/charade/zoomx/wiki/Manual>`__ for more informations

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