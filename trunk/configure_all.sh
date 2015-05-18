#!/bin/bash

INSTDIR=${HOME}/local

pushd lib_aux

#unpack gsl
tar xvzf gsl-1.16.tar.gz
#compile and install gsl
pushd gsl-1.16
./configure --prefix=${INSTDIR}
if [ $? -eq 0 ]
then
    make
    if [ $? -eq 0 ]
    then
	make install
    else
	echo "Make Failed"
	exit 2
    fi
else
    echo "Configure Failed"
    exit 1
fi
popd
#delete gsl
rm -rf gsl-1.16

tar xvzf cfitsio.tar.gz
#compile and install cfitiso
pushd cfitsio
./configure --prefix=${INSTDIR}
if [ $? -eq 0 ]
then
    make
    if [ $? -eq 0 ]
    then
        make install
    else
        echo "Make Failed"
        exit 2
    fi
else
    echo "Configure Failed"
    exit 1
fi
popd
#delete cfitsio
rm -rf cfitsio

tar xvzf CCfits-2.4.tar.gz
#compile and install CCfits                                                                                                                                                                             
pushd CCfits
./configure --prefix=${INSTDIR} --with-cfitsio=${INSTDIR}
if [ $? -eq 0 ]
then
    make
    if [ $? -eq 0 ]
    then
        make install
    else
        echo "Make Failed"
        exit 2
    fi
else
    echo "Configure Failed"
    exit 1
fi
popd
#delete CCfits
rm -rf CCfits

tar xvzf pyfits-3.3.tar.gz
#compile and install pyfits
pushd pyfits-3.3
python setup.py install --root=${INSTDIR}
if [ $? -ne 0 ]
then
    echo "Make Failed"
    exit 1
fi
popd
#delete pyfits
rm -rf pyfits-3.3

#exit lib_aux
popd

./configure --prefix=${INSTDIR} --LDFLAGS=-L${INSTDIR}/lib CPPFLAGS=-I${INSTDIR}/include
if [ $?-eq 0 ]
then
    make
    if [ $? -eq0 ]
    then
        make install
    else
	echo "Make Failed"
	exit 2
    fi
else
    echo "Configure Failed"
    exit 1
fi

echo "Compilation and Installationg Successful"
echo "For Full functionality add"
echo '  export SURVEYSIMPATH='${INSTDIR}
echo '  export PYTHONPATH=${PYTHONPATH}:'${INSTDIR}'/opt/TWWfsw/python27/lib/python2.7/site-packages:'${INSTDIR}'/surveysim/python'
echo "To your initialization file, or make sure to manually set them before execution." 
echo "These should only be necessary if you do not have root access, otherwise all default"
echo "install paths should be used"