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

Staged germline/somatic analysis with zoomxRun
--------------------------------------------------

$ZOOMXPKG/scripts/zoomxRun

  .. code:: bash

    #!/bin/bash
    
    #this script pipelines ZoomX germline and somatic call analysis
    
    #germline test configuration as in germline3.cfg
    #bam from 10X Chromium data of Chr22:22M-24M of NA12878
    
    #somatic test configuration as in somatic7.cfg
    #bam from 10X Chromium data of Chr8:98M-100M of a matched MetB7175 and Norm7176 sample
    
    #tested on BioLinux v8 Virtual Machine on running on AWS
    #hardware configuration: cpu=1,mem=2GB,disk=25GB
    #set ZOOMXPKG to your downloaded or cloned package path
    
    if [[ -z $1 ]]; then
      echo "#usage:"
      echo "#zoomxRun <stage> $cfgfile"  
      echo "#to auto generate a list of command lines, do following"
      echo "export cfgfile=some/cfgfile"
      echo "zoomxRun"
      echo "#------------------------------------------"
      echo "#individual/germline sample analysis stages"
      echo "#------------------------------------------"
      echo "zoomxRun caseCov $cfgfile       #prepare coverage tracks"  
      echo "zoomxRun caseParse $cfgfile     #parse Bam from long-ranger"
      echo "zoomxRun caseScan $cfgfile      #do grid Poisson field scan"
      echo "zoomxRun caseDistRef $cfgfile   #refine large-scale (>200K) intra-chrosomsome i.e. dist events"  
      echo "zoomxRun caseTransRef $cfgfile  #refine large-scale inter-chrosomsome i.e. trans events"
      echo "zoomxRun caseCisRef $cfgfile    #refine mid-scale (50K-200K) intra-chrosomsome i.e. cis events"  
      echo "zoomxRun caseDistPlot $cfgfile  #plot dist events in case (optional)"  
      echo "zoomxRun caseTransPlot $cfgfile #plot trans events in case (optional)"  
      echo "zoomxRun caseCisPlot $cfgfile   #plot cis events in case (optional)"  
      echo "#------------------------------------------"
      echo "#matched samples analysis stages           "
      echo "#------------------------------------------"
      echo "zoomxRun controlCov $cfgfile"
      echo "zoomxRun controlParse $cfgfile"  
      echo "zoomxRun controlScan $cfgfile"  
      echo "zoomxRun controlCisRef $cfgfile"  
      echo "zoomxRun controlDistRef $cfgfile"  
      echo "zoomxRun controlTransRef $cfgfile"
      echo "zoomxRun controlDistPlot $cfgfile  #plot dist events in control (optional)"  
      echo "zoomxRun controlTransPlot $cfgfile #plot trans events in control (optional)"  
      echo "zoomxRun controlCisPlot $cfgfile   #plot cis events in control (optional)"  
      echo "zoomxRun distDiffMix $cfgfile    #find somatic dist events assuming matched is mixture"  
      echo "zoomxRun transDiffMix $cfgfile   #find somatic trans events assuming matched is mixture"  
      echo "zoomxRun cisDiffMix $cfgfile     #find somatic cis events assuming matched is mixture"  
      echo "zoomxRun distDiffPure $cfgfile   #find somatic dist events assuming matched is pure"  
      echo "zoomxRun transDiffPure $cfgfile  #find somatic trans events assuming matched is pure"  
      echo "zoomxRun cisDiffPure $cfgfile    #find somatic cis events assuming matched is pure"  
      echo "zoomxRun distShare $cfgfile      #find shared dist events between the two matched"  
      echo "zoomxRun transShare $cfgfile     #find shared trans events between the two matched"  
      echo "zoomxRun cisShare $cfgfile       #find shared cis events between the two matched"  
      echo "#add keyword dry to any <stage> to preview the auto-generated commands" 
      exit
    fi
    
    testSubstring () { # if $1 is a substring of $2
      echo $2 | grep -E -o $1
    }
    
    ZoomXRunOp=$1
    ZoomxCfgFile=$2
    dryMode=$( testSubstring "dry" $ZoomXRunOp )
    addHeader='sed "1ichr1\tstart1\tend1\tchr2\tstart2\tend2\tmolcnt\trdpIM\trdpFF\trdpRR"'
    source $ZoomxCfgFile
    
    echo "{"
    echo "  'ZOOMXPATH': '$ZOOMXPATH',"
    echo "  'Task': 'ZoomX Calling'," 
    echo "  'Note': 'to run: change configurations in \$ZoomxCfgFile and specify \$ZoomXRunOp accordingly',"
    echo "  'mydate': '$mydate',"
    echo "  'ZoomXRunOp': '$ZoomXRunOp',"
    echo "  'ZoomxCfgFile': '$ZoomxCfgFile',"
    echo "  'dryMode': '$dryMode',"
    echo "  'ref': '$ref',"
    echo "  'grid': '$grid',"
    echo "  'covGrid': '$covGrid',"
    echo "  'outputPf': '$outputPf',"
    echo "  'distSize': $distSize,"
    echo "  'cisSize': $cisSize,"
    echo "  'controlPf': '$controlPf',"
    echo "  'controlBam': '$controlBam',"
    echo "  'controlExp': $controlExp,"
    echo "  'controlMax': $controlMax,"
    echo "  'casePf': '$casePf',"
    echo "  'caseBam': '$caseBam',"
    echo "  'caseExp': $caseExp,"
    echo "  'caseMax': $caseMax,"
    echo "  'maxPlot': $maxPlot,"
    
    #process control files;
    if [[ ! -z $( testSubstring "control" $ZoomXRunOp ) ]]; then # define common control variables;
      log=$outputPf/$controlBam.log
      bed=$outputPf/$controlBam.multiFragBar.bed.a0
      clean=$bed.clean
      filter=$clean.bed
      h5=$outputPf/$controlBam.multiFragBar.h5.a0
      trans=$outputPf/$controlBam.jct.bedpe.tra   #output of inter chromosome junctions
      distal=$outputPf/$controlBam.jct.bedpe.dis  #output of distal (>200 KB) junctions
      cis=$outputPf/$controlBam.jct.bedpe.cis     #output of intra  (>50 KB) junctions
      cfg=$outputPf/$controlBam.gridScan.cfg.json        #output of gridScan configuraiton
      echo "  'controlOption': {"
      echo "    'log': '$log', 'bed': '$bed', 'clean':'$clean', 'filter':'$filter', 'h5': '$h5', 'trans': '$trans', 'distal': '$distal', 'cis':'$cis', 'cfg':'$cfg' },"
      echo "  'controlCmd': {"
    fi
    
    if [[ ! -z $( testSubstring "controlCov" $ZoomXRunOp ) ]]; then
      if [[ $caseCovFile == "none" ]]; then
        controlCovCMD="Error: need to specify a controlCovFile other than none in $ZoomxCfgFile!"
      else
        controlCovCMD="(time samtools bedcov $ref.bed.w100s100 $controlPf/$controlBam >$controlCovFile) 2>$log.$ZoomXRunOp.$mydate"
      fi
      if [[ -z $dryMode ]]; then echo $controlCovCMD | bash; fi
      echo "    controlCovCMD': '$controlCovCMD',"
    fi
    
    if [[ ! -z $( testSubstring "controlParse" $ZoomXRunOp ) ]]; then
      parse10xCMD="PYTHONPATH=$ZOOMXPATH  python -m parse10x -o $outputPf/$controlBam -l 2 $controlPf/$controlBam >$log.$ZoomXRunOp.$mydate 2>&1"        #run parse10x module
      if [[ -z $dryMode ]]; then echo $parse10xCMD | bash; fi
      echo "    parse10xCMD': '$parse10xCMD',"
    fi
    
    if [[ ! -z $( testSubstring "controlScan" $ZoomXRunOp ) ]]; then
      if [[ $controlCovFile != "none" && -e $controlCovFile ]]; then  #covFile is specified and generated
        filterCMD="$ZOOMXPATH/../scripts/filterMolecule.R $controlCovFile $bed >$log.$ZoomXRunOp.$mydate 2>&1" # apply filtering
      else #covFile is not specified or generated (missing?)
        filterCMD="(cp $bed $clean && touch $filter) >$log.$ZoomXRunOp.$mydate 2>&1" #create dummy filter files
      fi
      if [[ -z $dryMode ]]; then echo $filterCMD | bash; fi
      echo "    filterCMD': '$filterCMD',"
      gridScanCMD="PYTHONPATH=$ZOOMXPATH  python -m gridScan -n 3 $grid $clean $h5 $outputPf/$controlBam >>$log.$ZoomXRunOp.$mydate 2>&1" #run gridScan module on Chr22 only (50 Mbp)
      if [[ -z $dryMode ]]; then echo $gridScanCMD | bash; fi
      echo "    'gridScanCMD': '$gridScanCMD'," 
      cisJctCMD="awk '{if(\$1 == \$4 && \$5 - \$2>$cisSize && \$5 - \$2<=$distSize) print}' $outputPf/$controlBam.jct.bedpe >$cis 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $cisJctCMD | bash; fi
      echo "    'cisJctCMD':'$cisJctCMD'," 
      distJctCMD="awk '{if(\$1 == \$4 && \$5 - \$2>$distSize) print}' $outputPf/$controlBam.jct.bedpe >$distal 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $distJctCMD | bash; fi
      echo "    'distJctCMD': '$distJctCMD',"
      transJctCMD="awk '{if(\$1 != \$4 && \$1 < \$4) print}' $outputPf/$controlBam.jct.bedpe >$trans 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $transJctCMD | bash; fi
      echo "    'transJctCMD': '$transJctCMD',"
      IMrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$controlBam.IM.bedpe >$outputPf/$controlBam.geq$cisSize.IM.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $IMrdpCMD | bash; fi
      echo "    'IMrdpCMD':'$IMrdpCMD',"
      FFrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$controlBam.FF.bedpe >$outputPf/$controlBam.geq$cisSize.FF.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $FFrdpCMD | bash; fi
      echo "    'FFrdpCMD':'$FFrdpCMD',"
      RRrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$controlBam.RR.bedpe >$outputPf/$controlBam.geq$cisSize.RR.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $RRrdpCMD | bash; fi
      echo "    'RRrdpCMD':'$RRrdpCMD',"
    fi
      
    if [[ ! -z $( testSubstring "controlCisRef" $ZoomXRunOp ) ]]; then
      cisRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $controlMax -e $controlExp -f $controlExp -o $outputPf/$controlBam.cis.E$controlExp -i $outputPf/$controlBam.geq$cisSize.IM.bedpe -j $outputPf/$controlBam.geq$cisSize.FF.bedpe -k $outputPf/$controlBam.geq$cisSize.RR.bedpe $controlPf/$controlBam $cis $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $cisRefCMD | bash; fi
      echo "    'cisRefCMD':'$cisRefCMD',"
    fi
      
    if [[ ! -z $( testSubstring "controlDistRef" $ZoomXRunOp ) ]]; then
      distRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $controlMax -e $controlExp -f $controlExp -o $outputPf/$controlBam.dist.E$controlExp -i $outputPf/$controlBam.geq$cisSize.IM.bedpe -j $outputPf/$controlBam.geq$cisSize.FF.bedpe -k $outputPf/$controlBam.geq$cisSize.RR.bedpe $controlPf/$controlBam $distal $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $distRefCMD | bash; fi
      echo "    'distRefCMD':'$distRefCMD',"
    fi
      
    if [[ ! -z $( testSubstring "controlTransRef" $ZoomXRunOp ) ]]; then
      transRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $controlMax -e $controlExp -f $controlExp -o $outputPf/$controlBam.tra.E$controlExp -i $outputPf/$controlBam.geq$cisSize.IM.bedpe -j $outputPf/$controlBam.geq$cisSize.FF.bedpe -k $outputPf/$controlBam.geq$cisSize.RR.bedpe $controlPf/$controlBam $trans $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $transRefCMD | bash; fi
      echo "    'transRefCMD':'$transRefCMD',"
    fi
    
    #process case files
    if [[ ! -z $( testSubstring "case" $ZoomXRunOp ) ]]; then # define common case variables;
      log=$outputPf/$caseBam.log
      bed=$outputPf/$caseBam.multiFragBar.bed.a0
      clean=$bed.clean
      filter=$clean.bed
      h5=$outputPf/$caseBam.multiFragBar.h5.a0
      trans=$outputPf/$caseBam.jct.bedpe.tra   #output of inter chromosome junctions
      distal=$outputPf/$caseBam.jct.bedpe.dis  #output of distal (>200 KB) junctions
      cis=$outputPf/$caseBam.jct.bedpe.cis     #output of intra  (>50 KB) junctions
      cfg=$outputPf/$caseBam.gridScan.cfg.json        #output of gridScan configuraiton
      echo "  'caseOption': {"
      echo "    'log': '$log', 'bed': '$bed', 'clean':'$clean', 'filter':'$filter', 'h5': '$h5', 'trans': '$trans', 'distal': '$distal', 'cis':'$cis', 'cfg':'$cfg' },"
      echo "  'caseCmd': {"
    fi
    
    if [[ ! -z $( testSubstring "caseCov" $ZoomXRunOp ) ]]; then
      if [[ $caseCovFile == "none" ]]; then
        caseCovCMD="Error: need to specify a caseCovFile other than none in $ZoomxCfgFile!"
      else
        caseCovCMD="(time samtools bedcov $ref.bed.w100s100 $casePf/$caseBam >$caseCovFile) 2>$log.$ZoomXRunOp.$mydate"
      fi
      if [[ -z $dryMode ]]; then echo $caseCovCMD | bash; fi
      echo "    caseCovCMD': '$caseCovCMD',"
    fi
    
    if [[ ! -z $( testSubstring "caseParse" $ZoomXRunOp ) ]]; then
      parse10xCMD="PYTHONPATH=$ZOOMXPATH  python -m parse10x -c $caseCovFile -o $outputPf/$caseBam -l 2 $casePf/$caseBam >$log.$ZoomXRunOp.$mydate 2>&1"        #run parse10x module
      if [[ -z $dryMode ]]; then echo $parse10xCMD | bash; fi
      echo "    parse10xCMD': '$parse10xCMD',"
    fi
    
    if [[ ! -z $( testSubstring "caseScan" $ZoomXRunOp ) ]]; then
      if [[ $caseCovFile != "none" && -e $caseCovFile ]]; then  #covFile is specified and generated
        filterCMD="$ZOOMXPATH/../scripts/filterMolecule.R $caseCovFile $bed >$log.$ZoomXRunOp.$mydate 2>&1" # apply filtering
      else #covFile is not specified or generated (missing?)
        filterCMD="cp $bed $clean && touch $filter >$log.$ZoomXRunOp.$mydate 2>&1" #create dummy filter files
      fi
      if [[ -z $dryMode ]]; then echo $filterCMD | bash; fi
      echo "    filterCMD': '$filterCMD',"
      gridScanCMD="PYTHONPATH=$ZOOMXPATH  python -m gridScan -n 3 $grid $clean $h5 $outputPf/$caseBam >>$log.$ZoomXRunOp.$mydate 2>&1" #run gridScan module on Chr22 only (50 Mbp)
      if [[ -z $dryMode ]]; then echo $gridScanCMD | bash; fi
      echo "    'gridScanCMD': '$gridScanCMD'," 
      cisJctCMD="awk '{if(\$1 == \$4 && \$5 - \$2>$cisSize && \$5 - \$2<=$distSize) print}' $outputPf/$caseBam.jct.bedpe >$cis 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $cisJctCMD | bash; fi
      echo "    'cisJctCMD':'$cisJctCMD'," 
      distJctCMD="awk '{if(\$1 == \$4 && \$5 - \$2>$distSize) print}' $outputPf/$caseBam.jct.bedpe >$distal 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $distJctCMD | bash; fi
      echo "    'distJctCMD': '$distJctCMD',"
      transJctCMD="awk '{if(\$1 != \$4 && \$1 < \$4) print}' $outputPf/$caseBam.jct.bedpe >$trans 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $transJctCMD | bash; fi
      echo "    'transJctCMD': '$transJctCMD',"
      IMrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$caseBam.IM.bedpe >$outputPf/$caseBam.geq$cisSize.IM.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $IMrdpCMD | bash; fi
      echo "    'IMrdpCMD':'$IMrdpCMD',"
      FFrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$caseBam.FF.bedpe >$outputPf/$caseBam.geq$cisSize.FF.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $FFrdpCMD | bash; fi
      echo "    'FFrdpCMD':'$FFrdpCMD',"
      RRrdpCMD="awk '{if( \$4 != \$1 || (\$5-\$2)>$cisSize || (\$5-\$2)<-$cisSize) print}' $outputPf/$caseBam.RR.bedpe >$outputPf/$caseBam.geq$cisSize.RR.bedpe 2>>$log.$ZoomXRunOp.$mydate"
      if [[ -z $dryMode ]]; then echo $RRrdpCMD | bash; fi
      echo "    'RRrdpCMD':'$RRrdpCMD',"
    fi
      
    if [[ ! -z $( testSubstring "caseCisRef" $ZoomXRunOp ) ]]; then
      cisRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $caseMax -e $caseExp -f $caseExp -o $outputPf/$caseBam.cis.E$caseExp -i $outputPf/$caseBam.geq$cisSize.IM.bedpe -j $outputPf/$caseBam.geq$cisSize.FF.bedpe -k $outputPf/$caseBam.geq$cisSize.RR.bedpe $casePf/$caseBam $cis $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $cisRefCMD | bash; fi
      echo "    'cisRefCMD':'$cisRefCMD',"
    fi
     
    if [[ ! -z $( testSubstring "caseDistRef" $ZoomXRunOp ) ]]; then
      distRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $caseMax -e $caseExp -f $caseExp -o $outputPf/$caseBam.dist.E$caseExp -i $outputPf/$caseBam.geq$cisSize.IM.bedpe -j $outputPf/$caseBam.geq$cisSize.FF.bedpe -k $outputPf/$caseBam.geq$cisSize.RR.bedpe $casePf/$caseBam $distal $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $distRefCMD | bash; fi
      echo "    'distRefCMD':'$distRefCMD',"
    fi
      
    if [[ ! -z $( testSubstring "caseTransRef" $ZoomXRunOp ) ]]; then
      transRefCMD="PYTHONPATH=$ZOOMXPATH  python -m refineBreak -n 8 -m $caseMax -e $caseExp -f $caseExp -o $outputPf/$caseBam.tra.E$caseExp -i $outputPf/$caseBam.geq$cisSize.IM.bedpe -j $outputPf/$caseBam.geq$cisSize.FF.bedpe -k $outputPf/$caseBam.geq$cisSize.RR.bedpe $casePf/$caseBam $trans $cfg >$log.$ZoomXRunOp.$mydate 2>&1"
      if [[ -z $dryMode ]]; then echo $transRefCMD | bash; fi
      echo "    'transRefCMD':'$transRefCMD',"
    fi
    
    #plotting of Diff or Share Molecules, moleFile is predefined in previous Diff and Share stages
    if [[ ! -z $(testSubstring "Plot" $ZoomXRunOp) ]]; then 
      if [[ ! -z $(testSubstring "case" $ZoomXRunOp) && ! -z $(testSubstring "Dist" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$caseBam.dist.E$caseExp.refined" )
      fi 
      if [[ ! -z $(testSubstring "case" $ZoomXRunOp) && ! -z $(testSubstring "Trans" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$caseBam.trans.E$caseExp.refined" )
      fi 
      if [[ ! -z $(testSubstring "case" $ZoomXRunOp) && ! -z $(testSubstring "Cis" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$caseBam.cis.E$caseExp.refined" )
      fi 
      if [[ ! -z $(testSubstring "control" $ZoomXRunOp) && ! -z $(testSubstring "Dist" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$controlBam.dist.E$controlExp.refined" )
      fi 
      if [[ ! -z $(testSubstring "control" $ZoomXRunOp) && ! -z $(testSubstring "Trans" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$controlBam.trans.E$controlExp.refined" )
      fi 
      if [[ ! -z $(testSubstring "control" $ZoomXRunOp) && ! -z $(testSubstring "Cis" $ZoomXRunOp) ]]; then
        moleFile=$( echo "$outputPf/$controlBam.cis.E$controlExp.refined" )
      fi 
      echo "  'getAndPlotCMD':["
      if [[ -e $moleFile ]]; then
        log="$moleFile.log.$ZoomXRunOp.$mydate"
        maxSize=$(cat $moleFile | wc -l)
        maxSize=$(expr $(($maxPlot>$maxSize?$maxSize:$maxPlot)) - 1)
        for i in `seq 1 $maxSize`; do
          getMolCMD="$ZOOMXPATH/../scripts/getMolecule.sh $moleFile $outputPf/$caseBam $i $moleFile.Plot >>$log 2>&1"
          if [[ -z $dryMode ]]; then echo $getMolCMD | bash; fi
          plotMolCMD="$ZOOMXPATH/../scripts/plotMolecule.R $isUCSC $moleFile.Plot_$i $outputPf/$caseBam.gridScan.cfg.json >>$log 2>&1" 
          if [[ -z $dryMode ]]; then echo $plotMolCMD | bash; fi
          echo "    [ '$getMolCMD', '$plotMolCMD' ],"
        done
      fi
      echo "  ]"
    fi
    
    if [[ ! -z $( testSubstring "case" $ZoomXRunOp ) || ! -z $( testSubstring "control" $ZoomXRunOp ) ]]; then
      echo "  }" #closing caseCmd or controlCmd
    fi
    
    #do somatic differentiation, contamination mode, assuming normal tissus is a mix too
    if [[ ! -z $( testSubstring "cisDiffMix" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.cis.E$caseExp"_"$controlBam.cis.E$controlExp.somatic )
      echo "    'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      cisDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.cis.E$caseExp.refined -b $outputPf/$controlBam.cis.E$controlExp.refined >$outputPf/$caseBam.cis.E$caseExp\_$controlBam.cis.E$controlExp.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $cisDiffCMD | bash; fi
      echo "  'cisDiffCMD':'$cisDiffCMD',"
    fi
    
    if [[ ! -z $( testSubstring "distDiffMix" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.dist.E$caseExp"_"$controlBam.dist.E$controlExp.somatic )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      distDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.dist.E$caseExp.refined -b $outputPf/$controlBam.dist.E$controlExp.refined >$outputPf/$caseBam.dist.E$caseExp\_$controlBam.dist.E$controlExp.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $distDiffCMD | bash; fi
      echo "  'distDiffCMD':'$distDiffCMD',"
    fi
    
    if [[ ! -z $( testSubstring "transDiffMix" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.tra.E$caseExp"_"$controlBam.tra.E$controlExp.somatic )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      transDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.tra.E$caseExp.refined -b $outputPf/$controlBam.tra.E$controlExp.refined >$outputPf/$caseBam.tra.E$caseExp\_$controlBam.tra.E$controlExp.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $transDiffCMD | bash; fi
      echo "  'transDiffCMD':'$transDiffCMD',"
    fi
    
    #do somatic differentiation, pure mode, assuming normal tissus is a pure germline
    if [[ ! -z $( testSubstring "cisDiffPure" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.cis.E$caseExp"_"$controlBam.cis.E0.somatic )
      echo "    'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      cisDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.cis.E$caseExp.refined -b $outputPf/$controlBam.jct.bedpe.cis >$outputPf/$caseBam.cis.E$caseExp\_$controlBam.cis.E0.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $cisDiffCMD | bash; fi
      echo "  'cisDiffCMD':'$cisDiffCMD',"
    fi
    
    if [[ ! -z $( testSubstring "distDiffPure" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.dist.E$caseExp"_"$controlBam.dist.E0.somatic )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      distDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.dist.E$caseExp.refined -b $outputPf/$controlBam.jct.bedpe.dis >$outputPf/$caseBam.dist.E$caseExp\_$controlBam.dist.E0.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $distDiffCMD | bash; fi
      echo "  'distDiffCMD':'$distDiffCMD',"
    fi
    
    if [[ ! -z $( testSubstring "transDiffPure" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.tra.E$caseExp"_"$controlBam.tra.E0.somatic )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      transDiffCMD="bedtools pairtopair -is -type notboth -a $outputPf/$caseBam.tra.E$caseExp.refined -b $outputPf/$controlBam.jct.bedpe.tra >$outputPf/$caseBam.tra.E$caseExp\_$controlBam.tra.E0.somatic 2>$log"
      if [[ -z $dryMode ]]; then echo $transDiffCMD | bash; fi
      echo "  'transDiffCMD':'$transDiffCMD',"
    fi
    
    #do germline sharing
    if [[ ! -z $( testSubstring "cisShare" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.cis.E$caseExp"_"$controlBam.cis.E$controlExp.germline )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      cisShareCMD="bedtools pairtopair -is -type both -a $outputPf/$caseBam.cis.E$caseExp.refined -b $outputPf/$controlBam.cis.E$controlExp.refined | cut -f1,2,3,4,5,6,7,8,9,10 | $addHeader >$outputPf/$caseBam.cis.E$caseExp\_$controlBam.cis.E$controlExp.germline 2>$log"
      if [[ -z $dryMode ]]; then echo $cisShareCMD | bash; fi
      echo "  'cisShareCMD':'$cisShareCMD',"
    fi
    
    if [[ ! -z $( testSubstring "distShare" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.dist.E$caseExp"_"$controlBam.dist.E$controlExp.germline )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      distShareCMD="bedtools pairtopair -is -type both -a $outputPf/$caseBam.dist.E$caseExp.refined -b $outputPf/$controlBam.dist.E$controlExp.refined | cut -f1,2,3,4,5,6,7,8,9,10 | $addHeader >$outputPf/$caseBam.dist.E$caseExp\_$controlBam.dist.E$controlExp.germline 2>$log"
      if [[ -z $dryMode ]]; then echo $distShareCMD | bash; fi
      echo "  'distShareCMD':'$distShareCMD',"
    fi
    
    if [[ ! -z $( testSubstring "transShare" $ZoomXRunOp ) ]]; then
      moleFile=$( echo $outputPf/$caseBam.tra.E$caseExp"_"$controlBam.tra.E$controlExp.germline )
      echo "  'moleFile':'$moleFile',"
      log="$moleFile.log.$ZoomXRunOp.$mydate"
      transShareCMD="bedtools pairtopair -is -type both -a $outputPf/$caseBam.tra.E$caseExp.refined -b $outputPf/$controlBam.tra.E$controlExp.refined | cut -f1,2,3,4,5,6,7,8,9,10 | $addHeader >$outputPf/$caseBam.tra.E$caseExp\_$controlBam.tra.E$controlExp.germline 2>$log"
      if [[ -z $dryMode ]]; then echo $transShareCMD | bash; fi
      echo "  'transShareCMD':'$transShareCMD',"
    fi
    
    #plotting of Diff or Share Molecules, moleFile is predefined in previous Diff and Share stages
    if [[ ! -z $(testSubstring "Diff" $ZoomXRunOp) || ! -z $(testSubstring "Share" $ZoomXRunOp) ]]; then 
      echo "  'getAndPlotCMD':["
      if [[ -e $moleFile ]]; then
        log="$moleFile.log.$ZoomXRunOp.$mydate"
        maxSize=$(cat $moleFile | wc -l)
        maxSize=$(expr $(($maxPlot>$maxSize?$maxSize:$maxPlot)) - 1)
        for i in `seq 1 $maxSize`; do
          getMolCMD="$ZOOMXPATH/../scripts/getMolecule.sh $moleFile $outputPf/$caseBam $i $moleFile.Plot >>$log 2>&1"
          if [[ -z $dryMode ]]; then echo $getMolCMD | bash; fi
          plotMolCMD="$ZOOMXPATH/../scripts/plotMolecule.R $isUCSC $moleFile.Plot_$i $outputPf/$caseBam.gridScan.cfg.json >>$log 2>&1" 
          if [[ -z $dryMode ]]; then echo $plotMolCMD | bash; fi
          echo "    [ '$getMolCMD', '$plotMolCMD' ],"
        done
      fi
      echo "  ]"
    fi
    
     #here we iterate through following steps for each junction file, handling max maxPlot events
     #  to get plotMolecule.R to work will need R>3 and R packages optparse, rjson, ggbio, ggplot2, GenomicRanges, RColorBrewer
     #  the packages should be automatically installed at first run of plotMolecule.R
     #  if not, you can do it mannually
    
    echo "}" #closing global
