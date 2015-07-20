#!/bin/bash

pushd obs
./update_makefile.sh
popd

pushd templates
./update_makefile.sh
popd

pushd model
./update_makefile.sh
popd
