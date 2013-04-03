#!/bin/sh
# Download and build a package from GitHub

PROJ=spagedi
PROJ_DISTS=SPAGeDi-1*
GH=https://github.com/reedacartwright/spagedi/archive/
MAKE=make
CMAKE=cmake

build_mingw32=
build_m32=
build_archive=HEAD
for build_option; do
	case $build_option in
	--mingw32)
		build_mingw32=yes ;;
	--m32)
		build_m32=yes ;;
	*)
		build_archive=$build_option ;;
	esac
done

RELENG_DIR=`mktemp -q -d -t ${PROJ}-releng.XXX`

echo Using $RELENG_DIR as a temporary build directory ...

if [ $? -ne 0 ]; then
        echo "$0: Can't create temp directory, exiting..."
        exit 1
fi

DEST_DIR=`pwd`
cd $RELENG_DIR || exit 1

archive=${PROJ}.tar.gz

echo 
echo Downloading SPAGeDi Archive from GitHub ...

wget -nv -O ${archive} "${GH}${build_archive}.tar.gz" || exit 1

echo 
echo Identifying Archive ...

version=`tar -tf ${archive} --exclude='*/*' | sed "s/^${PROJ}-\(.*\)\/\$/\1/"`
version_len=`expr $version : '.*'`

if [ $version_len -eq 40 ]; then
	version=`echo $version | cut -c1-7`
	echo "    Archive" is commit $version.
else
	echo "    Archive" is tag $version.
fi

echo 
echo Extracting Archive ...

tar xzf ${archive} && cd ${PROJ}-* || exit 1
mkdir -p build && cd build

echo
echo Configuring Source ...

if test $build_mingw32; then
	$CMAKE .. -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_TOOLCHAIN_FILE="../releng/mingw32.cmake" \
		-DUSE_STATIC_LIBS=on \
		-DSPAGEDI_VERSION="${version}"
elif test $build_m32; then
	$CMAKE .. -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_C_FLAGS=-m32 \
		-DUSE_STATIC_LIBS=on \
		-DSPAGEDI_VERSION="${version}"
else
	$CMAKE .. -DCMAKE_BUILD_TYPE=Release \
	-DUSE_STATIC_LIBS=on \
	-DSPAGEDI_VERSION="${version}"
fi

echo
echo Building Binary Package ...
$MAKE new_package || exit 1

echo
echo Build Source Package ...
$MAKE new_package_source || exit 1

echo
echo Moving distribution packages ...

mv -v ${PROJ_DISTS} ${DEST_DIR}

echo
echo Cleaning up ...

cd $DEST_DIR
rm -rf $RELENG_DIR

