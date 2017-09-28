PHASE="mus.txt"
SCEL=0
CONF=0

for MU in $(cat $PHASE); do
	NEWDIR=mu$MU
	mkdir $NEWDIR
	cd $NEWDIR
	shadowdir ../../..
    rm $(cat ../removable.txt)
	cd ..
	NEWSETTINGS="$NEWDIR/monte_settings.json"
	sed "s/RRRR/$MU/" monte_settings.json > $NEWSETTINGS
	sed -i "s/SSSS/$SCEL/" $NEWSETTINGS
	sed -i "s/CCCC/$CONF/" $NEWSETTINGS
	cp static.pbs $NEWDIR
done
