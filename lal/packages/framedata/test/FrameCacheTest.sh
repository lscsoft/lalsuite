#!/bin/sh

./FrameCacheTest || exit
sed '/^[X|L]/d' catalog.test | cmp - catalog.out
exit
