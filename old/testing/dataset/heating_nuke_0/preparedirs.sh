#!/bin/bash

PHASE="phase3.txt"
SCEL=49
CONF=4

for MU in $(cat $PHASE); do
	NEWDIR=mu$MU
	mkdir $NEWDIR
	cd $NEWDIR
	shadowdir ../../..
	cd ..
	NEWSETTINGS="$NEWDIR/monte_settings.json"
	sed "s/RRRR/$MU/" monte_settings.json > $NEWSETTINGS
	sed -i "s/SSSS/$SCEL/" $NEWSETTINGS
	sed -i "s/CCCC/$CONF/" $NEWSETTINGS
	cp static.pbs $NEWDIR
done
