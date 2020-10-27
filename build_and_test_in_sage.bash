#!/bin/bash

set -e

cp -r /multipolynomial-bases .
cd multipolynomial-bases

echo
echo
echo "==================================================================="
echo "Install make"
sudo apt-get update && sudo apt-get install make


echo
echo
echo "==================================================================="
echo "Insall git"
sudo apt-get install -y git

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

echo
echo
echo "==================================================================="
echo "Deploying the doc"
set -e
ABS_DEPLOY_KEY="`pwd`/.travis_ci_gh_pages_deploy_key"
if [[ -r "$ABS_DEPLOY_KEY" ]]; then
    echo "Deployment key exists, attempting to upload"
    DEPLOY_DOC_TO_REPOSITORY="VivianePons/multipolynomial-bases"
    chmod 600 "$ABS_DEPLOY_KEY"
    export GIT_SSH_COMMAND="ssh -v -i $ABS_DEPLOY_KEY -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
    rm -Rf gh-pages
    echo "Attempting git clone, repo $DEPLOY_DOC_TO_REPOSITORY"
    git clone --depth 1 git@github.com:${DEPLOY_DOC_TO_REPOSITORY}.git --depth 1 --branch=gh-pages gh-pages
    BUILT_DOCS_DIR=`cd /home/sage/multipolynomial-bases/docs/build/html && pwd`
    echo "cd gh-pages"
    cd gh-pages
    echo "rm -Rf ./${DEPLOY_DOC_TO_DIRECTORY}/*"
    rm -Rf ./${DEPLOY_DOC_TO_DIRECTORY}/*
    echo "mkdir -p ./${DEPLOY_DOC_TO_DIRECTORY}"
    mkdir -p ./${DEPLOY_DOC_TO_DIRECTORY}
    echo "cp -R $BUILT_DOCS_DIR/* ./${DEPLOY_DOC_TO_DIRECTORY}/"
    cp -R $BUILT_DOCS_DIR/* ./${DEPLOY_DOC_TO_DIRECTORY}/
    git add --all .
    git config user.name "Travis CI"
    git config user.email "nobody@example.org"
    if git commit -m "Automatic upload of documentation built from ${TRAVIS_COMMIT}"; then
    echo "Attempting git push"
	git push origin gh-pages
    fi
fi

