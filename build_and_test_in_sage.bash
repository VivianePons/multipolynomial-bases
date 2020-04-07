#!/bin/bash

set -e

cp -r /multipolynomial-bases .
cd multipolynomial-bases

RUN apt install make

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

