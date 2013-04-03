#!/bin/sh
# Download and build a package from GitHub

CMAKE=cmake

# Process command line arguments
TEMP=`getopt -o t: --long toolchain: -n $0 -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
build_toolchain=
while true ; do
	case "$1" in
		-t|--toolchain) build_toolchain=$2 ; shift 2 ;;
		--) shift ; break ;;
		*) echo "Internal error!" ; exit 1 ;;
	esac
done

# Determine the archive tag
build_archive=${1-HEAD}
build_args="-DRELENG_TAG=${build_archive}"

# if toolchain is m32 or m64 use the flags and not a toolchain
if [ -n "${build_toolchain}" ]; then
	case "${build_toolchain}" in
		m32|M32) build_args="${build_args} -DRELENG_M32=on" ;;
		m64|M64) build_args="${build_args} -DRELENG_M64=on" ;;
		*)       build_args="${build_args} -DRELENG_TOOLCHAIN=${build_toolchain}" ;;
	esac
fi

$CMAKE $build_args -P releng.cmake

