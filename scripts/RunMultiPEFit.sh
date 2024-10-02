#!/bin/bash
echo "Checking directory: $1"
if [ "${1: -4}" == "root"  ]
	then
		nameroot=${1}
		echo "Wave data exists"
		if [ -e ${nameroot} ]
			then
				echo "Converted data exists."
				echo ${nameroot}
				root -l -b -q 'MultiPEFit.C("'${nameroot}'","'${2}'")'
		fi
	else
		echo "Wave data does not exist, skipping directory."
fi