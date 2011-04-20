#!/bin/sh

PROJ=spagedi
PROJ_DISTS=spagedi-1*
MAKE=make
CMAKE=cmake
REPOS=`svn info | grep URL: | perl -pe 's!^URL: (.+)/releng$!$1!'`

build_mingw32=
build_m32=
for build_option; do
	case $build_option in
	--mingw32)
		build_mingw32=yes ;;
	--m32)
		build_m32=yes ;;
	esac
done

RELENG_DIR=`mktemp -q -d -t ${PROJ}-releng.XXX`
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

if test -n $build_mingw32; then
	$CMAKE $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_TOOLCHAIN_FILE="${SOURCE_DIR}/releng/i586-mingw32msvc.cmake" \
		-DUSE_STATIC_LIBS=on
elif test -n $build_m32; then
	$CMAKE $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release
else
	$CMAKE $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release
fi
$MAKE
$MAKE package
$MAKE package_source

echo
echo Moving distribution packages ...

mv -v SPAGeDi-1* $DEST_DIR

echo
echo Cleaning up ...

cd $DEST_DIR
rm -rf $RELENG_DIR

