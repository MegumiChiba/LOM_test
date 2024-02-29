#!/bin/bash
#echo "wb_num: $1"
#mask=$((2**$1))
echo "pwren_mask: $mask"
export LOMTEST="$HOME/LOMtest/STM32Workspace_1220/src/tools/python/LOMTest"

python3 $HOME/LOMtest/fh_server/scripts/pmt_hv_enable.py -w 0
#wblist=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
wblist=(15)
#vlist=(76 78 80 82 84 86 88 90 92 94 96)
vlist=(80 82 84)
for nwb in "${wblist[@]}"; do
  #$HOME/LOMtest/LOMpre.sh
  python ~/tmp/mon_pressure.py
  echo "nwb: $nwb"
  mask=$((2**$nwb))
  echo "mask: $mask"
  #python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config spestatus.cfg --ofile lom_setup.txt --nobatch
  for vc in "${vlist[@]}"; do
      $HOME/LOMtest/LOMpre.sh
      echo "Vc: $vc"
      echo "$nwb $mask"
      python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config ./config/gain$vc.cfg --ntosend_total 10000000 --ntosend_batch 10 --runtime 60 --ofile $HOME/LOMtest/lom-test-data/gaintest/lom$nwb.gain$vc.dat --binary #--noboot
      FILE="$HOME/LOMtest/lom-test-data/gaintest/lom$nwb.gain$vc.dat"
      if [ -s ${FILE} ]; then
          echo "datasize OK" 
      else
         $HOME/LOMtest/LOMpre.sh 
          echo "Retry data taking"
          python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask $mask --wb_num $nwb --wub_config ./config/gain$vc.cfg --ntosend_total 10000000 --ntosend_batch 10 --runtime 60 --ofile $HOME/LOMtest/lom-test-data/gaintest/lom$nwb.gain$vc.dat --binary --noboot
      fi
  done

#python3 $HOME/LOMtest/lom-test-data/plot-binary-data.py lom_spe_$1 lom_spe_$1.dat

  python3 $LOMTEST/lom_single_base_batch_data.py --devmode --loglevel info --host=localhost --port=5000 --pwren_mask 0 --wub_config spe.cfg --ofile lom_base_off.txt --nobatch --nosetup
  python plot-binary-gain.py /home/icehap/LOMtest/lom-test-data/gaintest/lom$nwb.gain $nwb
  python read_gain_pdf.py /home/icehap/LOMtest/lom-test-data/gaintest/lom$nwb.gain /home/icehap/LOMtest/lom-test-data/gaintest/pdfs/lom$nwb.gain $nwb
done
