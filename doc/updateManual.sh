#!/bin/bash
export ZOOMXPKG=$HOME/setup/zoomx                #where the zoomx package is located
export ZOOMXPATH=$ZOOMXPKG/zoomx                 #where the python libs of zoomx are located
export SCRIPTPATH=$ZOOMXPKG/scripts              #where the uitlity scripts zoomx are located
export PYTHONPATH=$ZOOMXPATH
export PATH=$SCRIPTPATH:$PATH
cat Manual.master.rst | awk '/REPLACE_WITH_ZOOMX_PARSE10X_HELP/{system("python -m parse10x -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_GRIDSCAN_HELP/{system("python -m gridScan -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_REFINEBREAK_HELP/{system("python -m refineBreak -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_FILTERMOLECULE_HELP/{system("filterMolecule.R -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_PLOTMOLECULE_HELP/{system("plotMolecule.R -h | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_ZOOMXRUN_HELP/{system("zoomxRun | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_ZOOMXSTAGE_HELP/{system("zoomxStage | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_GETMOLECULE_HELP/{system("getMolecule.sh | sed -e \"s/^/    /g\"");next}1' | awk '/REPLACE_WITH_ZOOMX_BEDPE2LOUPE_HELP/{system("bedpe2loupe.sh | sed -e \"s/^/    /g\"");next}1' >Manual.rst
