#!/bin/bash

set -e

cp -r /multipolynomial-bases .
cd multipolynomial-bases

echo
echo
echo "==================================================================="
echo "Install make"
sudo apt-get update && sudo apt-get install make


# echo
# echo
# echo "==================================================================="
# echo "Insall git"
# sudo apt-get install -y git

echo
echo
echo "==================================================================="
echo "SageMath version"
sage --version


echo
echo
echo "==================================================================="
echo "Installing Multipolynomial bases"
make install

echo
echo
echo "==================================================================="
echo "Testing Myltupolynomial bases"
make test


echo
echo
echo "==================================================================="
echo "Building the doc"
make doc


