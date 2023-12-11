#!/bin/bash

while true; do
    MEMK=`cat /proc/meminfo | grep MemAvailable | rev | cut -d" " -f2 | rev`
    MEMG=`expr ${MEMK} \* 1024 | numfmt --to=iec`
    SWAPK=`vmstat -s | grep "used swap" | awk '{$1=$1};1' | cut -d" " -f1`
    SWAPG=`expr ${SWAPK} \* 1024 | numfmt --to=iec`
    DISKROOT=`df --output=avail -h / | awk '{$1=$1};1' | grep -v Avail`
    DISKDATA=`df --output=avail -h /home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/flybase_metta | awk '{$1=$1};1' | grep -v Avail`
    echo "RAM: ${MEMG} | Swap: ${SWAPG} | Root: ${DISKROOT} | Data: ${DISKDATA}"
    sleep 60
done
