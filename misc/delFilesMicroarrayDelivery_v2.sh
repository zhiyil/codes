#!/bin/bash
set -euo pipefail
 
tree -fi | egrep "_WF\.xls|[Cc]y3\.txt$|\.apg$|\.tab$|Single.+\.xls$|Serv.+\.xlsx?$|anl$|Report\.xls$|txt$|MultiArray\\\ Analysis\.xls$" > delList
echo -e "files to be deleted:\n"
cat delList
sed -i '' 's/\\//g' delList

IFS=$'\n'
while true; do
    read -p "ok to delete these files? y or n" yn
    case $yn in
        [Yy]* ) for file in `cat delList`; do rm ${file}; done; rm delList; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no...";;
    esac
done

