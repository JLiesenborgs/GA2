#!/bin/bash -e

if [ -z "$NUM" ] ; then
	NUM=20
fi

FUNCTIONNAMES=" \
f1 \
f2 \
f3 \
f5 \
f6 \
f7 \
f8 \
f9_k4 \
f9_k8 \
f11_30 \
f11_100 \
f12_10 \
f12_30 \
f13_20 \
f13_100 \
f14_20 \
f14_100 \
f15_30 \
f15_100 \
f16 \
f17 \
f18 \
f19_0.5 \
f19_1 \
f20 \
f21_2 \
f21_3 \
f21_4 \
f22_5 \
f22_8 \
f22_10 \
f23_2 \
f23_3 \
f23_4 \
f24_5 \
f24_6 \
f24_7 \
f25 \
f26 \
f27 \
f28_1 \
f28_2 \
f28_3 \
f28_4 \
f28_5 \
f28_6 \
f29 \
"

if [ $# -lt 1 ] ; then
	echo "Specify the location of the debasictest executable"
	exit -1
fi

EXE="$1"
shift

if [ $# -gt 0 ] ; then
	FUNCTIONNAMES="$1"
fi

for FUNCTIONNAME in $FUNCTIONNAMES ; do

	echo
	echo "$FUNCTIONNAME =========================================================="
	for i in `seq 1 $NUM` ; do 
		"$EXE" "$FUNCTIONNAME" ; 
	done | tee tmplog
	num=`grep evaluations tmplog | wc -l`

	FAILED=`grep -v "eval:" tmplog |wc -l`
	FAILEDSTR=""
	B0=""
	B1=""
	if [ "$FAILED" -gt 0 ] ; then
		B0="("
		B1=")"
		FAILEDSTR="Failed: $FAILED/$NUM"
	fi

	NFE=$( ( echo -n "(0" ; for i in `grep "eval:" tmplog | cut -f 12 -d " "` ; do echo -n "+$i" ; done ; echo ")/$NUM" ) | bc -l )
	NFE=`printf "%.02f" $NFE`

	echo "$FUNCTIONNAME Average evaluations: ${B0}${NFE}${B1} $FAILEDSTR"
done
