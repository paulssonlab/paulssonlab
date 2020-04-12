#!/bin/sh
bold="$(tput bold)"
normal="$(tput sgr0)"
echo "🐑  ${bold}Cloning git@github.com:paulssonlab/paulssonlab.git into ${normal}${PWD}/paulssonlab"
git clone git@github.com:paulssonlab/paulssonlab.git
status=$?
if [ $status != 0 ] ; then
    echo "❌  ${bold}Could not clone repo, aborting.${normal}"
    exit 1
fi
cd paulssonlab
bin/init-repo
status=$?
if [ $status != 0 ] ; then
    echo "❌  ${bold}Could not initialize repo. Please fix issue and rerun bin/init-repo.${normal}"
    exit 1
fi
