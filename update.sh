#!/bin/bash

git pull
pushd Source
./configure
make
sudo make install
make clean
popd

echo "SurveySim updated and built"