.. contents:: Navigator

Examples
========

Download and unzip ZoomX package or clone the source repository to $ZOOMXPKG.
ZoomX are most versatile if used from Python and R scripts directly.
To easy the initial learning, we also include a BASH script zoomxRun.
The script stages common analysis tasks with easily editable configure files. 

An example ZoomX germline analysis
--------------------------------------------------
To run this example, run the 'germline3.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz>`__

$ZOOMXPKG/test/germline3.cfg

  .. code:: bash

    #!/bin/bash
    #if zoomx python entry-points (parse10x, gridScan, refineBreak) 
    #and utility scripts are not in your PATH
    #ZOOMXPKG has to be set here 
    ZOOMXPKG=$HOME/setup/zoomx                #where the zoomx package is located
    ZOOMXPATH=$ZOOMXPKG/zoomx                 #where the python libs of zoomx are located
    SCRIPTPATH=$ZOOMXPKG/scripts              #where the uitlity scripts zoomx are located
    PATH=$SCRIPTPATH:$PATH
    mydate=$(LC_ALL=en_US.utf8 date +%Y%b%d)
    ref=human_g1k_v37.ucsc.fasta #prefix to grid file, actual fasta file not needed
    grid=$ref.bed.w10Ks5K
    covGrid=$ref.bed.w100s100
    outputPf=.
    distSize=200000
    cisSize=50000
    casePf=.
    caseBam=NA12878.wg12.chr22-22M-24M.bam
    caseExp=45
    caseMax=1000
    caseCovFile=$casePf/$caseBam.cov100
    isUCSC="--isUCSC" #turn on --isUCSC if chrX not X
    maxPlot=100

$ZOOMXPKG/test/germline3.run

  .. code:: bash

    #!/bin/bash
    if [[ -z $1 ]]; then
      cfgfile="germline3.cfg"
    else
      cfgfile=$1
    fi
    
    source $cfgfile &&
    zoomxRun caseCov $cfgfile >$cfgfile.caseCov.json.$mydate &&
    echo "$cfgfile.caseCov.json.$mydate" &&
    zoomxRun caseParse $cfgfile >$cfgfile.caseParse.json.$mydate &&
    echo "$cfgfile.caseParse.json.$mydate" &&
    zoomxRun caseScan $cfgfile >$cfgfile.caseScan.json.$mydate &&
    echo "$cfgfile.caseScan.json.$mydate" &&
    zoomxRun caseDistRef $cfgfile >$cfgfile.caseDistRef.json.$mydate &&
    echo "$cfgfile.caseDistRef.json.$mydate" &&
    zoomxRun caseTransRef $cfgfile >$cfgfile.caseTransRef.json.$mydate &&
    echo "$cfgfile.caseTransRef.json.$mydate" &&
    zoomxRun caseCisRef $cfgfile >$cfgfile.caseCisRef.json.$mydate &&
    echo "$cfgfile.caseCisRef.json.$mydate" &&
    zoomxRun caseDistPlot $cfgfile >$cfgfile.caseDistPlot.json.$mydate &&
    echo "$cfgfile.caseDistPlot.json.$mydate" &&
    zoomxRun caseTransPlot $cfgfile >$cfgfile.caseTransPlot.json.$mydate &&
    echo "$cfgfile.caseTransPlot.json.$mydate" &&
    zoomxRun caseCisPlot $cfgfile >$cfgfile.caseCisPlot.json.$mydate &&
    echo "$cfgfile.caseCisPlot.json.$mydate"

An example ZoomX somatic analysis
--------------------------------------------------

To run this example, run the 'somatic7.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz>`__

test/somatic7.cfg

  .. code:: bash

    #!/bin/bash
    #if zoomx python entry-points (parse10x, gridScan, refineBreak)
    #and utility scripts are not in your PATH
    #ZOOMXPKG has to be set here
    ZOOMXPKG=$HOME/setup/zoomx                #where the zoomx package is located
    ZOOMXPATH=$ZOOMXPKG/zoomx                 #where the python libs of zoomx are located
    SCRIPTPATH=$ZOOMXPKG/scripts              #where the uitlity scripts zoomx are located
    PATH=$SCRIPTPATH:$PATH
    mydate=$(LC_ALL=en_US.utf8 date +%Y%b%d)
    ref=human_g1k_v37.fasta                   #prefix to grid file, actual fasta file not needed
    grid=$ref.bed.w10Ks5K                     #grid file
    covGrid=$ref.bed.w100s100                 #grid file for coverage
    outputPf=.
    distSize=200000
    cisSize=50000
    runOption=""
    controlPf=.
    controlBam=Norm7176.wg12.chr8-98M-100M.bam
    controlExp=10
    controlMax=1000
    controlCovFile=$controlPf/$controlBam.cov100
    casePf=.
    caseBam=MetB7175.wg12.chr8-98M-100M.bam
    caseExp=25
    caseMax=100
    caseCovFile=$casePf/$caseBam.cov100
    isUCSC="" #turn on --isUCSC if chrX not X
    maxPlot=100

test/somatic7.run

  .. code:: bash

    #!/bin/bash
    if [[ -z $1 ]]; then
      cfgfile="somatic7.cfg"
    else
      cfgfile=$1
    fi
    
    source $cfgfile &&
    zoomxRun controlCov $cfgfile >$cfgfile.controlCov.json.$mydate &&
    echo $cfgfile.controlCov.json.$mydate &&
    zoomxRun caseCov $cfgfile >$cfgfile.caseCov.json.$mydate &&
    echo $cfgfile.caseCov.json.$mydate &&
    zoomxRun controlParse $cfgfile >$cfgfile.controlParse.json.$mydate &&
    echo $cfgfile.controlParse.json.$mydate &&
    zoomxRun caseParse $cfgfile >$cfgfile.caseParse.json.$mydate &&
    echo $cfgfile.caseParse.json.$mydate &&
    zoomxRun controlScan $cfgfile >$cfgfile.controlScan.json.$mydate &&
    echo $cfgfile.controlScan.json.$mydate &&
    zoomxRun caseScan $cfgfile >$cfgfile.caseScan.json.$mydate &&
    echo $cfgfile.caseScan.json.$mydate &&
    zoomxRun controlCisRef $cfgfile >$cfgfile.controlCisRef.json.$mydate &&
    echo $cfgfile.controlCisRef.json.$mydate &&
    zoomxRun controlDistRef $cfgfile >$cfgfile.controlDistRef.json.$mydate &&
    echo $cfgfile.controlDistRef.json.$mydate &&
    zoomxRun controlTransRef $cfgfile >$cfgfile.controlTransRef.json.$mydate &&
    echo $cfgfile.controlTransRef.json.$mydate &&
    zoomxRun caseCisRef $cfgfile >$cfgfile.caseCisRef.json.$mydate &&
    echo $cfgfile.caseCisRef.json.$mydate &&
    zoomxRun caseDistRef $cfgfile >$cfgfile.caseDistRef.json.$mydate &&
    echo $cfgfile.caseDistRef.json.$mydate &&
    zoomxRun caseTransRef $cfgfile >$cfgfile.caseTransRef.json.$mydate &&
    echo $cfgfile.caseTransRef.json.$mydate &&
    zoomxRun distDiffMix $cfgfile >$cfgfile.distDiffMix.json.$mydate &&
    echo $cfgfile.distDiffMix.json.$mydate &&
    zoomxRun cisDiffMix $cfgfile >$cfgfile.cisDiffMix.json.$mydate &&
    echo $cfgfile.cisDiffMix.json.$mydate &&
    zoomxRun transDiffMix $cfgfile >$cfgfile.transDiffMix.json.$mydate &&
    echo $cfgfile.transDiffMix.json.$mydate &&
    zoomxRun distDiffPure $cfgfile >$cfgfile.distDiffPure.json.$mydate &&
    echo $cfgfile.distDiffPure.json.$mydate &&
    zoomxRun cisDiffPure $cfgfile >$cfgfile.cisDiffPure.json.$mydate &&
    echo $cfgfile.cisDiffPure.json.$mydate &&
    zoomxRun transDiffPure $cfgfile >$cfgfile.transDiffPure.json.$mydate &&
    echo $cfgfile.transDiffPure.json.$mydate &&
    zoomxRun distShare $cfgfile >$cfgfile.distShare.json.$mydate &&
    echo $cfgfile.distShare.json.$mydate &&
    zoomxRun cisShare $cfgfile >$cfgfile.cisShare.json.$mydate &&
    echo $cfgfile.cisShare.json.$mydate &&
    zoomxRun transShare $cfgfile >$cfgfile.transShare.json.$mydate &&
    echo $cfgfile.transShare.json.$mydate
