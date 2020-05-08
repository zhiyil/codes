#!/bin/bash
 
cd "Standard Analysis"
echo -e "files to be deleted:\n"
find . -name "*WF.xls*" -exec ls {} \;
find . -name "*cy3.txt|*Cy3.txt" -exec ls {} \;
find . -name "*apg" -exec ls {} \;
find . -name "*.tab" -exec ls {} \;
find . -name "Serv*xls*" -exec ls {} \;
find . -name "Single*xls" -exec ls {} \;
 
while true; do
    read -p "ok to delete these files? y or n" yn
    case $yn in
        [Yy]* ) find . -name "*WF.xls*" -exec rm {} \;; find . -name "*cy3.txt|*Cy3.txt" -exec rm {} \;; find . -name "*apg" -exec rm {} \;; find . -name "*.tab" -exec rm {} \;; find . -name "Serv*xls*" -exec rm {} \;; find . -name "Single*xls" -exec rm {} \;; echo "Done!"; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no...";;
    esac
done
 
cd ../"In_Depth Analysis"
 
echo -e "files to be deleted:\n"
find . -name "*anl" -exec ls {} \;
find . -name "*Report.xls" -exec ls {} \;
find . -name "*txt" -exec ls {} \;
find . -name "*MultiArray*xls" -exec ls {} \;
 
while true; do
    read -p "ok to delete these files? y or n" yn
    case $yn in
        [Yy]* ) find . -name "*anl" -exec rm {} \;; find . -name "*Report.xls" -exec rm {} \;; find . -name "*txt" -exec rm {} \;; find . -name "*MultiArray*xls" -exec rm {} \;; echo "Done!"; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no...";;
    esac
done