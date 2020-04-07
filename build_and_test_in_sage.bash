#!/bin/bash

set -e

cp -r /multipolynomial-bases .
cd multipolynomial-bases

echo
echo
echo "==================================================================="
echo "SageMath version"
sage --version


echo
echo
echo "==================================================================="
echo "Installing Multipolynomial basis"
make install

echo
echo
echo "==================================================================="
echo "Testing SnapPy"
make test

