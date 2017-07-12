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

echo    # move to a new line
read -p "Do you want to add the gin bins and libs to your environment paths in .bashrc? [Y/n]" -n 1 -r
echo    # move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

cat <<- EOF  >> ~/.bashrc

# added by gin installer
export CPLUS_INCLUDE_PATH=$INSTALL_PATH/include:\${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=$INSTALL_PATH/lib:\${LIBRARY_PATH}
export PATH=$INSTALL_PATH/bin:\${PATH}
EOF

source ~/.bashrc

fi