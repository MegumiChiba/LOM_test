#!/bin/bash
export LOMTEST="$HOME/LOMtest/STM32Workspace_1220/src/tools/python/LOMTest"
python3 $HOME/LOMtest/fh_server/scripts/pmt_hv_enable.py -w 0
wblist=(11 12) #0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
index=(0 3 5 7) # 5 6 7 8 9)

for i in "${index[@]}"; do
    echo "Loop $i reset"
    #$HOME/LOMtest/LOMpre.sh
    for nwb in "${wblist[@]}"; do
        $HOME/LOMtest/LOMpre.sh
        echo "Loop: $i"
        echo "nwb: $nwb"
        mask=$((2**$nwb))
        echo "mask: $mask"
        python ~/tmp/mon_pressure.py
        #python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config ./config/spestatus.cfg --ofile lom_setup.txt --nobatch
        python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config ./config/dark$nwb.cfg --ntosend_total 20000000 --ntosend_batch 30 --runtime 300 --ofile $HOME/LOMtest/lom-test-data/darktest/lom$nwb.dark$i.dat --binary #--noboot
        FILE="$HOME/LOMtest/lom-test-data/darktest/lom$nwb.dark$i.dat"
        if [ -s ${FILE} ]; then
            echo "datasize OK" 
        else 
            echo "Retry data taking"
            python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config ./config/dark$nwb.cfg --ntosend_total 20000000 --ntosend_batch 10 --runtime 300 --ofile $HOME/LOMtest/lom-test-data/darktest/lom$nwb.dark$i.dat --binary #--noboot
        fi

    #done
    python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask 0 --wub_config spe.cfg --ofile lom_base_off.txt --nobatch --nosetup
    python plot-binary-dark.py $HOME/LOMtest/lom-test-data/darktest/pdfs/log/lom$nwb.dark$i $nwb $HOME/LOMtest/lom-test-data/darktest/lom$nwb.dark$i
    done
done
for nwb in "${wblist[@]}"; do
    python read_dark_pdf.py $HOME/LOMtest/lom-test-data/darktest/lom$nwb.dark $HOME/LOMtest/lom-test-data/darktest/pdfs/lom$nwb.dark $nwb
done
