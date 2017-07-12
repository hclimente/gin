#!env bash

cd gin
cmake -DCMAKE_INSTALL_PREFIX:PATH=build .
make all install

echo    # move to a new line
read -p "Do you want to add the gin bins and libs to your environment paths in .bashrc? [Y/n]" -n 1 -r
echo    # move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

cat <<- EOF  >> ~/.bashrc

# added by gin installer
export CPLUS_INCLUDE_PATH=`pwd`/build/include:\${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=`pwd`/build/lib:\${LIBRARY_PATH}
export PATH=`pwd`/build/bin:\${PATH}
EOF

source ~/.bashrc

fi