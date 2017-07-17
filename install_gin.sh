#!/bin/bash

# crash in case of error
set -e

INSTALL_PATH=`pwd`/gin/build

if [ ! -z "$1" ]
  then
    INSTALL_PATH=$1
fi

cd gin
cmake -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_PATH .
make all install

export CPLUS_INCLUDE_PATH=$INSTALL_PATH/include:${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=$INSTALL_PATH/lib:${LIBRARY_PATH}
export PATH=$INSTALL_PATH/bin:${PATH}

cat <<- EOF

Copy to your .bashrc to add gin and the executables to your PATHs:

export CPLUS_INCLUDE_PATH=$INSTALL_PATH/include:\${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=$INSTALL_PATH/lib:\${LIBRARY_PATH}
export PATH=$INSTALL_PATH/bin:\${PATH}
EOF