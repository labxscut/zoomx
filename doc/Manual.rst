.. contents:: Navigator

Manual
========

1. Short Manuals for Python and R scripts are available with '-h' option
-------------------------------------------------------------------------

$ZOOMXPKG/zoomx/parse10x.py

.. code:: bash

    usage: parse10x.py [-h] [-c COVFILE] [-l LONGRANGERVERSION] [-f SETALPHA]
                       [-b MAXREAD] [-o OUTPREFIX] [-a]
                       bam
    
    positional arguments:
      bam                   sorted 10x bamfile
    
    optional arguments:
      -h, --help            show this help message and exit
      -c COVFILE, --covFile COVFILE
                            100bp window coverage file
      -l LONGRANGERVERSION, --longRangerVersion LONGRANGERVERSION
                            v>1.2 will use 10X's MI tag to infer molecules if
                            specified and available
      -f SETALPHA, --setAlpha SETALPHA
                            alpha i.e. max molecule extension basepair, override
                            by 0 in case of using MI tag
      -b MAXREAD, --maxRead MAXREAD
                            max number of read to process
      -o OUTPREFIX, --outPrefix OUTPREFIX
                            outputPrefix
      -a, --debug           in debug mode

$ZOOMXPKG/zoomx/gridScan.py

.. code:: bash

    usage: gridScan.py [-h] [-b GRIDFILE2] [-s SKIPFILE1] [-t SKIPFILE2]
                       [-e EXPECTFALSEDISCOVER] [-x MAXDUPLICATION]
                       [-y MAXJUNCTION] [-p ANEUPLOIDYALLOWANCE] [-n NUMSUB]
                       [-m MAXTAGFRACTION] [-g GENOMESIZE] [-r GRIDSIZE] [-a]
                       [-d TEMPDIR]
                       gridFile1 fragFile hdf5File outputPrefix
    
    positional arguments:
      gridFile1             1st set of grids BED formatted
      fragFile              10x fragments BED formatted
      hdf5File              10x statistics HDF5 formatted
      outputPrefix          output prefix
    
    optional arguments:
      -h, --help            show this help message and exit
      -b GRIDFILE2, --gridFile2 GRIDFILE2
                            2nd set of grids BED formatted (optional) (default:
                            same)
      -s SKIPFILE1, --skipFile1 SKIPFILE1
                            1st skipped grid regions BED formatted (optional)
                            (default: none)
      -t SKIPFILE2, --skipFile2 SKIPFILE2
                            2nd skipped grid regions BED formatted (optional)
                            (default: same)
      -e EXPECTFALSEDISCOVER, --expectFalseDiscover EXPECTFALSEDISCOVER
                            expected false discovery allowed (default: 10)
      -x MAXDUPLICATION, --maxDuplication MAXDUPLICATION
                            maximum expected number of >gridSize duplications
                            (default: 2147483647)
      -y MAXJUNCTION, --maxJunction MAXJUNCTION
                            maximum expected number of >gridSize junctions
                            (default: 2147483647)
      -p ANEUPLOIDYALLOWANCE, --aneuploidyAllowance ANEUPLOIDYALLOWANCE
                            adjust background rate by known aneuploidy (default:
                            4)
      -n NUMSUB, --numSub NUMSUB
                            number of submatrix used for huge sparse matrices dot
                            multiplication, memory=1/n^2 (default: 2)
      -m MAXTAGFRACTION, --maxTagFraction MAXTAGFRACTION
                            max tag fraction to be considerrred (default: -1)
      -g GENOMESIZE, --genomeSize GENOMESIZE
                            bp size of the whole genome (default: 3200000000)
      -r GRIDSIZE, --gridSize GRIDSIZE
                            bp size of the unit grid (default: 10000)
      -a, --debug           in debug mode (default: False)
      -d TEMPDIR, --tempdir TEMPDIR
                            root dir for keeping temporary files, default (last 4
                            digits random): /home/lixia/zoomx_jcXq

$ZOOMXPKG/zoomx/gridScan.py

.. code:: bash

    usage: refineBreak.py [-h] [-c COVFILE] [-m MAXNUM] [-e ALTEXPECT]
                          [-f FORCEDEXPECT] [-i IMBEDPEFILE] [-j FFBEDPEFILE]
                          [-k RRBEDPEFILE] [-o OUTPREFIX] [-a] [-b] [-n NUMPROCS]
                          [-d TEMPDIR]
                          bam call cfg
    
    positional arguments:
      bam                   sorted 10x bamfile, in Longranger BAM format
      call                  SV call bedpefile, in bedpe format
      cfg                   parameter config file, in JSON format
    
    optional arguments:
      -h, --help            show this help message and exit
      -c COVFILE, --covFile COVFILE
                            100bp window coverage file
      -m MAXNUM, --maxNum MAXNUM
                            max number of events to process (reduce this if time
                            constrained)
      -e ALTEXPECT, --altExpect ALTEXPECT
                            minimum mixture fraction for alternative, specify a
                            percentage between 10-100, default 25
      -f FORCEDEXPECT, --forcedExpect FORCEDEXPECT
                            forced minimum molecule for events
      -i IMBEDPEFILE, --imBedpeFile IMBEDPEFILE
                            impp FR read file in bedpe, default=bam.IM.bedpe
      -j FFBEDPEFILE, --ffBedpeFile FFBEDPEFILE
                            disc FF read file in bedpe, default=bam.FF.bedpe
      -k RRBEDPEFILE, --rrBedpeFile RRBEDPEFILE
                            splt RR read file in bedpe, default=bam.RR.bedpe
      -o OUTPREFIX, --outPrefix OUTPREFIX
                            outputPrefix
      -a, --debug           in debug mode
      -b, --mergeOnly       in mergeOnly mode
      -n NUMPROCS, --numprocs NUMPROCS
                            number of parallele process
      -d TEMPDIR, --tempdir TEMPDIR
                            root dir for keeping temporary files, default (last 4
                            digits random): /home/lixia/zoomx_FiEn

$ZOOMXPKG/scripts/filterMolecule.R

.. code:: bash

    Usage: /home/lixia/setup/zoomx/scripts/filterMolecule.R [options] covFile moleFile
    
    
    Options:
    	-o OUTPREFIX, --outPrefix=OUTPREFIX
    		none=$molFile.clean [default none]
    
    	-x XMADEXCLUDE, --xMadExclude=XMADEXCLUDE
    		if coverage is x mad away from median, exclude regional molecules [default 1]
    
    	-y YMEDIANINCLUDE, --yMedianInclude=YMEDIANINCLUDE
    		if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default 6]
    
    	-z ZMADINCLUDE, --zMadInclude=ZMADINCLUDE
    		if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default 6]
    
    	-w WINDOWSINCLUDE, --windowsInclude=WINDOWSINCLUDE
    		if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default 35]
    
    	-q, --noQuiet
    		show verbose, [default FALSE]
    
    	-a, --debug
    		save debug, [default FALSE]
    
    	-h, --help
    		Show this help message and exit
    
    

$ZOOMXPKG/scripts/plotMolecule.R

.. code:: bash

    Usage: /home/lixia/setup/zoomx/scripts/plotMolecule.R [options] plotDir statsJSON
    
    
    Options:
    	-f FLUSHSD, --flushSD=FLUSHSD
    		maximum of flushing sd [default 5000]
    
    	-p GAPSIZE, --gapSize=GAPSIZE
    		2 times the boundary padding bp [default 100000]
    
    	-s SCALESIZE, --scaleSize=SCALESIZE
    		boundary rounded to bp [default 10000]
    
    	-b BREAKCI, --breakCI=BREAKCI
    		limit of CI to estimate break point [default 0.95]
    
    	-u, --isUCSC
    		show verbose, [default FALSE]
    
    	-q, --noQuiet
    		show verbose, [default FALSE]
    
    	-a, --debug
    		save debug, [default FALSE]
    
    	-h, --help
    		Show this help message and exit
    
    

2. Short Manuals for Bash scripts are available with comment lines or call wihtout options  
-------------------------------------------------------------------------------------------

$ZOOMXPKG/scripts/zoomxRun

.. code:: bash

    #usage:
    #zoomxRun <stage> 
    #to auto generate a list of command lines, do following
    export cfgfile=some/cfgfile
    zoomxRun
    #------------------------------------------
    #individual/germline sample analysis stages
    #------------------------------------------
    zoomxRun caseCov        #prepare coverage tracks
    zoomxRun caseParse      #parse Bam from long-ranger
    zoomxRun caseScan       #do grid Poisson field scan
    zoomxRun caseDistRef    #refine large-scale (>200K) intra-chrosomsome i.e. dist events
    zoomxRun caseTransRef   #refine large-scale inter-chrosomsome i.e. trans events
    zoomxRun caseCisRef     #refine mid-scale (50K-200K) intra-chrosomsome i.e. cis events
    zoomxRun caseDistPlot   #plot dist events in case (optional)
    zoomxRun caseTransPlot  #plot trans events in case (optional)
    zoomxRun caseCisPlot    #plot cis events in case (optional)
    #------------------------------------------
    #matched samples analysis stages           
    #------------------------------------------
    zoomxRun controlCov 
    zoomxRun controlParse 
    zoomxRun controlScan 
    zoomxRun controlCisRef 
    zoomxRun controlDistRef 
    zoomxRun controlTransRef 
    zoomxRun controlDistPlot   #plot dist events in control (optional)
    zoomxRun controlTransPlot  #plot trans events in control (optional)
    zoomxRun controlCisPlot    #plot cis events in control (optional)
    zoomxRun distDiffMix     #find somatic dist events assuming matched is mixture
    zoomxRun transDiffMix    #find somatic trans events assuming matched is mixture
    zoomxRun cisDiffMix      #find somatic cis events assuming matched is mixture
    zoomxRun distDiffPure    #find somatic dist events assuming matched is pure
    zoomxRun transDiffPure   #find somatic trans events assuming matched is pure
    zoomxRun cisDiffPure     #find somatic cis events assuming matched is pure
    zoomxRun distShare       #find shared dist events between the two matched
    zoomxRun transShare      #find shared trans events between the two matched
    zoomxRun cisShare        #find shared cis events between the two matched
    #add keyword dry to any <stage> to preview the auto-generated commands

$ZOOMXPKG/zoomx/getMolecule.sh

.. code:: bash

    usage: getMolecule.sh inMoleFile inPrefix lineNum outDir
    example:
      getMolecule.sh NA12878.wg12.chr22.bam.dis.E40.refined NA12878.wg12.chr22.bam 1 testPlot
      #this gets singlue molecules supports 1st junction listed in NA12878.wg12.chr22.bam.dis.E40.refined
      #from NA12878.wg12.chr22.bam and store them in testPlot

$ZOOMXPKG/zoomx/bedpe2loupe.sh

.. code:: bash

    usage: bedpe2loupe.sh moleFile [chr]
    example: bedpe2loupe.sh NA12878.wg12.chr22.bam.dis.E40.refined
      convert bedpe coordinates to loupe coordinates. add [chr] to prefix if specified.

