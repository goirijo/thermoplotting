#!/bin/bash

for MUA in $(cat ./mua.txt); do
    for MUB in $(cat ./mub.txt); do
        NEWDIR=mu_${MUA}_${MUB}
        mkdir $NEWDIR
        sed "s/MUA/$MUA/" ./monte_settings.json > $NEWDIR/monte_settings.json
        sed -i "s/MUB/$MUB/" $NEWDIR/monte_settings.json
        sed "s/MUA/$MUA/" ./static.pbs > $NEWDIR/static.pbs
        sed -i "s/MUB/$MUB/" $NEWDIR/static.pbs
    done
done
