#!/bin/bash
echo "Measure pressure"
while :
do    
    python ~/tmp/mon_pressure.py
    python ~/tmp/offline_plot.py
    sleep 300
done
