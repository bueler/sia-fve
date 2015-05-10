#!/bin/bash

# generates a summary of history by looking in (e.g.) r900m?/history.txt and
# extracting some lines

set -e # exit on error
set -x

if [ -z ${FINE+x} ]; then
  RES=1800
  LEVELS="4 3 2 1"
else  # this case if FINE is *any* nonempty string
  RES=900
  LEVELS="5 4 3 2 1"
fi

# move old result out of way if present
OUT=study.${RES}m
touch $OUT
mv -f $OUT SAFE_$OUT

cathistory ()
{
head -n 1 $1 >> $2
cat $1 | grep "last successful value of eps" | sed 's/.* //g' >> $2
cat $1 | grep "maximum solution diffusivity" | sed 's/.* //g' >> $2
cat $1 | grep "total time" | sed 's/.* //g' >> $2
}

for LEV in $LEVELS; do
    cathistory r${RES}m${LEV}/history.txt $OUT
done

set +x

