#!/bin/bash

# Clone the remote dir to a new directory
# and create the standard directory structure

setupdir()
{
    echo "Creating dir $2"
    mkdir "$2"
    echo "Copying $src to $2"
    cp -r "$src"/* $2

    # Now create the required dirs
    mkdir "$2"/{build,logs,plots}
    mkdir -p "$2"/build/{diag,io,extern/unity,gfunc}
    mkdir -p "$2"/build/{ham_gen,params,tests,utils}

    # Symlink the data to the datadir
    # in the parent folder
    cd "$2"
    ln -s ../"$3" data
}

srcname="remote_clean"
src="../$srcname"
dirnamelist="flock leap zeal"
datadirname="data"

echo "Directory structure"
echo "|-- $srcname"
echo "|-- .."
echo "    |"
echo "    |-- datadir"
echo "    |"
echo "    |-- dirname"
echo "        |-- data -> ../datadir"
echo "        |-- build"
echo "        |-- src..."


# echo -n "Name of new dir: "
# read dirname
# echo -n "Name of data dir:"
# read datadirname
# echo ""
# setupdir $src $dirname $datadirname
# echo "Done"
# echo ""

for dir in $dirnamelist; do
    setupdir $src $dir $datadirname
    echo "Setup for $dir done"
    echo ""
done


exit 0