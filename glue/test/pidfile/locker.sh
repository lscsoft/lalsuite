#!/bin/sh

export PYTHONPATH=../..:$PYTHON_PATH

echo "=== Testing locking: this should succeed ===\n"

./locker.py

echo "\n=== Testing competing readers: this should produce a RuntimeError ===\n"

./locker.py 5 &
./locker.py
sleep 6

echo "\n=== Testing competing lock: this should produce a RuntimeError ===\n"

./locker.py 5 &
sleep 1
./locker.py
sleep 6

echo
