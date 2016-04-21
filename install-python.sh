#!/bin/sh
#
# script to install amuse build system from scratch
#
# author: Arjen van Elteren
# date  : 2009 - 05 -18
#
# Edit by TLRH to also install blas/lapack, pip, and use correct ssl version
# For BLAS/LAPACK: https://stackoverflow.com/questions/7496547

#APPVER=2.5.4
APPVER=2.7.9
#APPVER=2.6.5
#APPVER=2.7.1
APPFILE=Python-${APPVER}.tar.bz2
APPFILE=Python-${APPVER}.tgz
APP_DIR=Python-${APPVER}
URL=http://www.python.org/ftp/python/${APPVER}/${APPFILE}


# TODO: the installed version does not know its version?
# This throws warnings
OPENSSLVERSION="1.0.1r"
OPENSSLFILE=openssl-${OPENSSLVERSION}.tar.gz 
OPENSSLURL=http://mirrors.ibiblio.org/openssl/source/old/1.0.1/${OPENSSLFILE}
OPENSSLDIR=openssl-${OPENSSLVERSION}

BLASFILE=blas.tgz
BLASURL=http://www.netlib.org/blas/blas.tgz
LAPACKFILE=lapack.tgz
LAPACKURL=http://www.netlib.org/lapack/lapack.tgz

if [ -z ${PREFIX} ]; then
	echo The PREFIX variable is not set, please set it to an user directory
	exit 1
fi

if [ ! -d ${PREFIX} ]; then
	echo ${PREFIX} directory does not exists, please create it first!
	exit 1
 
fi

 

INSTALL_DIR=$PREFIX/install
mkdir $INSTALL_DIR
cd $INSTALL_DIR

DOWNLOAD_DIR=$INSTALL_DIR/_downloaded
BUILD_DIR=$INSTALL_DIR/_build
SOURCE_DIR=$INSTALL_DIR/_source

mkdir ${DOWNLOAD_DIR}
mkdir ${SOURCE_DIR}
mkdir ${BUILD_DIR}


echo "Downloading archive file..."
cd  ${DOWNLOAD_DIR}
if [ -e ${APPFILE} ] ; then
	echo "...File already downloaded";
else
	if which curl >/dev/null; then
		curl -L -O ${URL} ;
	else
		wget ${URL};		
	fi
fi

if [ -e ${OPENSSLFILE} ] ; then
	echo "...openssl file already downloaded";
else
	if which curl >/dev/null; then
		curl -L -O ${OPENSSLURL} ;
	else
		wget ${OPENSSLURL};
	fi
fi

if [ -e ${BLASFILE} ] ; then
	echo "...blas file already downloaded";
else
	if which curl >/dev/null; then
		curl -L -O ${BLASURL} ;
	else
		wget ${BLASURL};
	fi
fi

if [ -e ${LAPACKFILE} ] ; then
	echo "...lapack file already downloaded";
else
	if which curl >/dev/null; then
		curl -L -O ${LAPACKURL} ;
	else
		wget ${LAPACKURL};
	fi
fi

# TODO: add sqlite3 for jupyter notebook?
# https://stackoverflow.com/questions/8656158
# wget http://www.sqlite.org/sqlite-autoconf-3070900.tar.gz
# tar xvvf sqlite-autoconf-3070900.tar.gz
# cd sqlite-autoconf-3070900
# ./configure --prefix=~/applications
# make
# make install
# 
# cd ~/applications/src
# wget http://www.python.org/ftp/python/2.5.2/Python-2.5.2.tgz
# tar xvvf Python-2.5.2.tgz
# cd Python-2.5.2
# ./configure --prefix=~/applications
# make
# make install

cd ..
echo "Done"


cd ${SOURCE_DIR}
rm -Rf ${APP_DIR}
echo "Unpacking source files.."
tar -xf ${DOWNLOAD_DIR}/${APPFILE}
echo "..Done"

cd ${SOURCE_DIR}
rm -Rf ${OPENSSLDIR}
echo "Unpacking openssl source files.."
tar -xzf ${DOWNLOAD_DIR}/${OPENSSLFILE} || exit $?
echo "..Done"
cd ..

cd ${SOURCE_DIR}
rm -Rf ${BLASDIR}
echo "Unpacking blas source files.."
tar -xzf ${DOWNLOAD_DIR}/${BLASFILE} || exit $?
echo "..Done"
cd ..

cd ${SOURCE_DIR}
rm -Rf ${LAPACKDIR}
echo "Unpacking lapack source files.."
tar -xzf ${DOWNLOAD_DIR}/${LAPACKFILE} || exit $?
echo "..Done"
cd ..

# echo "Press enter to continue"
# read enterKey

echo "Building files.."

echo "Building openssl"

cd ${SOURCE_DIR}
cd ${OPENSSLDIR}


MACHINE=`(uname -m) 2>/dev/null`

./config \
    --prefix=${PREFIX}  \
    --openssldir=${PREFIX}/openssl \
    --shared

make -j8

make -j8 install

# echo "Press enter to continue"
# read enterKey

# blas/lapack
echo "Building BLAS"

cd ${SOURCE_DIR}
cd BLAS-*
gfortran -O3 -std=legacy -m64 -fno-second-underscore -fPIC -c *.f
ar r libfblas.a *.o
ranlib libfblas.a
rm -rf *.o
# export BLAS=~/src/BLAS-*/libfblas.a
mv libfblas.a ${PREFIX}/lib
export BLAS=${PREFIX}/lib/libfblas.a
ls -l ${BLAS}

# echo "Press enter to continue"
# read enterKey


echo "Building LAPACK"
cd ${SOURCE_DIR}
cd lapack-*
cp INSTALL/make.inc.gfortran make.inc  # On Linux with lapack-3.2.1 or newer
vi make.inc
make -j8 lapacklib
make -j8 clean
mv liblapack.a ${PREFIX}/lib
export LAPACK=${PREFIX}/lib/liblapack.a
ls -la ${LAPACK}

# echo "Press enter to continue"
# read enterKey

# python

cd ${BUILD_DIR}
rm -Rf ${APP_DIR}
mkdir ${APP_DIR}
cd ${APP_DIR}

UNAME=`uname`
if [ $UNAME == 'Darwin' ] ; then
	${SOURCE_DIR}/${APP_DIR}/configure \
		--with-dyld \
		--prefix=${PREFIX} \
	    --enable-unicode=ucs4\
		--program-suffix=.exe ;
#--enable-framework=${PREFIX}/Framework \
#--enable-universal-sdk \
else
	${SOURCE_DIR}/${APP_DIR}/configure \
		--prefix=${PREFIX} \
		--enable-shared \
	    --enable-unicode=ucs4\
		--program-suffix=.exe ;
fi

make -j8
make -j8 install
make -j8
make -j8 install
echo "..Done"

#. ~/.bash_functions
#which python
#setup_amuse
#which python
#exit 0
#echo "Install pip"
#cd ${SOURCE_DIR}
#wget https://bootstrap.pypa.io/get-pip.py
#python get-pip.py
echo "..Done"

#echo "Install scipy"
#pip install scipy
