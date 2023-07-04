#!/bin/bash

# Clone the remote dir to a new directory
# and create the standard directory structure

setupdir()
{
    echo "Creating dir $1"
    mkdir "$1"
    echo "Copying $src to $1"
    cp -r "$src"/* $1

    # Now create the required dirs
    mkdir "$1"/{build,logs,plots}
    mkdir -p "$1"/build/{diag,io,extern/unity,gfunc}
    mkdir -p "$1"/build/{ham_gen,params,tests,utils}

    # Symlink the data to the datadir
    # in the parent folder
    cd "$1"
    ln -s ../"$2" data
}

src="../remote_clean"

echo -n "Name of new dir: "
read dirname
echo "The data dir should be present"
echo "in the parent folder where"
echo "$dirname is to be created."
echo -n "Name of data dir:"
read datadirname

setupdir $dirname $datadirname
