#!/usr/bin/env bash

### Convert b/w image files (as e.g. a jpeg) to a .pgm file suitable to be used as geometry
# Note: forbidden geometries will NOT be checked/corrected! Also there is no guarantee
#   that rescaling will not introduce any forbidden geometry!
# Rescaling will be performed if optional width and height arguments are given
#
# Before usage you might need to make this file executable with:
#   chmod +x create-geometry.sh
#
# Usage:
#   ./create-geometry.sh FILE.jpg [width height]
#

USAGE="Usage:\n\t./create-geometry.sh FILE [width height]"
CONVERT_FLAGS="-posterize 2 -depth 1 -compress none -negate"

if [[ ! $1 ]]; then
    echo "ERROR: You need to pass a jpg file as argument!"
    echo -e ${USAGE}
    exit 1
fi

INFILE="$1"
BASENAME="${INFILE%.*}"
EXTENSION="${INFILE##*.}"
if [[ "${EXTENSION}" == "pgm" ]]; then
    echo "WARNING: a .pgm file was passed as input, do you want to perform the inverse convertion to .jpg?"
    echo -e "[y/N]:"
    read choice
    [[ "${choice}" != "y" ]] && exit 1
    OUTFILE="${BASENAME}.jpg"
    CONVERT_FLAGS="-negate"
else
    OUTFILE="${BASENAME}.pgm"
fi

if [[ $2 && $3 ]]; then
    WIDTH=$2
    HEIGHT=$3
    echo "INFO: rescaling the geometry to ${WIDTH}x${HEIGHT}"
    CONVERT_FLAGS="-scale ${WIDTH}x${HEIGHT} ${CONVERT_FLAGS}"
fi

echo "INFO: generating ${OUTFILE}"
convert ${CONVERT_FLAGS} ${INFILE} ${OUTFILE}; sed -i s/"255"/"1"/g ${OUTFILE}
exit $?

#eof
