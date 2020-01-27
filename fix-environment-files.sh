#!/bin/sh
retval=0
for file in "$@"; do
    if grep -q '^prefix:.*$' $file; then
        retval=1
        sed -i "" '/^prefix:.*$/d' $file
    fi
done
exit $retval
