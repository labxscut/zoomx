.. contents:: Navigator

Examples
========

ZoomX are most versatile if used from Python and R scripts directly.
To easy the initial learning, we also include a BASH script zoomxRun.
The script stages common analysis tasks with easily editable configure files. 

To run following examples:

  - Download and unzip ZoomX package or clone the source repository to $ZOOMXPKG.
  - Download and unzip the full test data set `zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__ and place the data into $ZOOMXPKG/test
  - Configure the $ZOOMXPKG/test/germline10.cfg and $ZOOMXPKG/test/somatic10.cfg files to reflect your $ZOOMXPKG
  - Change into $ZOOMXPKG/test and run the germline10.run and somatic10.run testing scripts
  
  .. code:: bash

    cd $ZOOMXPKG/test
    wget https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz -O zoomx.test.full.tgz
    tar -zxvf zoomx.test.full.tgz
    # edit $ZOOMXPKG/test/germline10.cfg  and $ZOOMXPKG/test/somatic10.cfg
    ./germline10.run
    ./somatic10.run

An example ZoomX germline analysis
--------------------------------------------------
To run this example, run the 'germline10.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$ZOOMXPKG/test/germline10.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINECFG

$ZOOMXPKG/test/germline10.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINERUN

An example ZoomX somatic analysis
--------------------------------------------------

To run this example, run the 'somatic10.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

test/somatic10.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICCFG

test/somatic10.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICRUN
