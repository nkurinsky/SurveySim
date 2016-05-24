#!/bin/bash

if [ $# -eq 1 ]
then
    echo "Making "$1" the current tag"
    rm current.tar.gz
    ln $1 current.tar.gz
else
    echo "Calling sequence: "$0" [current_tag]"
fi
