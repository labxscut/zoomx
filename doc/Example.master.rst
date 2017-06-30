.. contents:: Navigator

Examples
========

ZoomX are most versatile if used from Python and R scripts directly.
To easy the initial learning, we also include a BASH script zoomxRun.
The script stages common analysis tasks with easily editable configure files. 

To run following examples:

  - Download and unzip ZoomX package or clone the source repository to $ZOOMXPKG.
  - Download and unzip the full test data set `zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__ and place the data into $ZOOMXPKG/test
  - Configure the $ZOOMXPKG/test/germline9.cfg and $ZOOMXPKG/test/somatic9.cfg files to reflect your $ZOOMXPKG
  - Change into $ZOOMXPKG/test and run the germline9.run and somatic9.run testing scripts

An example ZoomX germline analysis
--------------------------------------------------
To run this example, run the 'germline9.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$ZOOMXPKG/test/germline9.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINECFG

$ZOOMXPKG/test/germline9.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINERUN

An example ZoomX somatic analysis
--------------------------------------------------

To run this example, run the 'somatic9.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

test/somatic9.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICCFG

test/somatic9.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICRUN
