#!/bin/sh

./FrameCacheTest || exit
sed '/^X/d' catalog.test | cmp - catalog.out
exit
