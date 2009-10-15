#!/bin/sh

./FrameCacheTest || exit
sed '/^[X|H|L]/d' catalog.test | cmp - catalog.out
exit
