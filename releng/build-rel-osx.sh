#!/bin/sh

MAKE=make
CMAKE=cmake
REPOS=svn://scit.us/documents/projects/spagedi

xcode=/Xcode2.5
xcodep=${xcode}/usr
PATH="${xcodep}/bin:${PATH}"; export PATH
CFLAGS="$CFLAGS -I${xcodep}/include"; export CFLAGS
CXXFLAGS="$CXXFLAGS -I${xcodep}/include"; export CXXFLAGS
LDFLAGS="-L${xcodep}/lib"; export LDFLAGS
CMAKE_OSX_SYSROOT="${xcode}/SDKs/MacOSX10.4u.sdk"; export CMAKE_OSX_SYSROOT

RELENG_DIR=`mktemp -q -d -t spagedi-releng`
if [ $? -ne 0 ]; then
        echo "$0: Can't create temp directory, exiting..."
        exit 1
fi

DEST_DIR=`pwd`
SOURCE_DIR="${RELENG_DIR}/source"
BUILD_DIR="${RELENG_DIR}/build"

svn co -q $REPOS $SOURCE_DIR || exit 1

mkdir $BUILD_DIR || exit 1
cd $BUILD_DIR || exit 1

$CMAKE $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=ppc;i386 \
  -DCPACK_SYSTEM_NAME=Darwin-universal \
  -DAPPLE_BUNDLE=ON
$MAKE
$MAKE package
$MAKE package_source

mv SPAGeDi-1* $DEST_DIR

cd $DEST_DIR
rm -rf $RELENG_DIR


