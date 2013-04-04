#!/bin/sh
# Download and build a package from GitHub

CMAKE=cmake

# Process command line arguments
build_toolchain=
while getopts t: name; do
	case $name in
		t) build_toolchain="$OPTARG" ;;
		?) printf "Usage: %s [-t toolchain] [tag]\n" $0
		   exit 2;;
	esac
done
shift `expr $OPTIND - 1`

# Determine the archive tag
build_archive=${1-HEAD}
build_args="-DRELENG_TAG=${build_archive}"

# if toolchain is m32 or m64 use the flags and not a toolchain
if [ ! -z "${build_toolchain}" ]; then
	case "${build_toolchain}" in
		m32|M32) build_args="${build_args} -DRELENG_M32=on" ;;
		m64|M64) build_args="${build_args} -DRELENG_M64=on" ;;
		*)       build_args="${build_args} -DRELENG_TOOLCHAIN=${build_toolchain}" ;;
	esac
fi

$CMAKE $build_args -P releng.cmake

