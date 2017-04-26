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
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$ZOOMXPKG/test/germline3.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINECFG

$ZOOMXPKG/test/germline3.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTGERMLINERUN

An example ZoomX somatic analysis
--------------------------------------------------

To run this example, run the 'somatic7.run' script in the 'test' subdir.
Full test data set is required, which can be downloaded at: 
`zoomx.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

test/somatic7.cfg

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICCFG

test/somatic7.run

  .. code:: bash

REPLACE_WITH_ZOOMXTESTSOMATICRUN

Staged germline/somatic analysis with zoomxRun
--------------------------------------------------

$ZOOMXPKG/scripts/zoomxRun

  .. code:: bash

REPLACE_WITH_ZOOMXRUN
