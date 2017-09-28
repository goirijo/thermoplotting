BACK=$(pwd)

for DIR in mu_*.*_*.*; do
    cd $DIR
    #sort by beta
    python $BACK/iteraverage.py | column -t | sort -n -k11 > tabulated_averages.txt
    cd $BACK
done
