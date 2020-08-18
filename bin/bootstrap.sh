#!/bin/sh
bold="$(tput bold)"
normal="$(tput sgr0)"
if [[ $HOSTNAME =~ login[0-9][0-9]\.o2\.rc\.hms\.harvard\.edu && ! $SLURMD_NODENAME =~ compute-.* ]]; then
    echo "⏱  ${bold}You are trying to run a long-running command on an O2 login node.${normal}"
    echo "You should probably abort and run"
    echo "    irun1"
    echo "to open an interactive session on a compute node. Then try again."
    read -n 1 -p "❌  ${bold}Do you want to abort? ${normal}[Y/n] " answer
    echo
    case ${answer:0:1} in
        n|N )
            exit 0
        ;;
        * )
            exit 1
        ;;
    esac
fi
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
