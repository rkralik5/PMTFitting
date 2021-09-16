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
				root -l 'Fit.C("'${nameroot}'")'
		fi
	else
		echo "Wave data does not exist, skipping directory."
fi

