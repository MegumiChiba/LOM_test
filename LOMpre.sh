#!/bin/bash
export PATH="$HOME/Workspace/DeggDaq/STF/stf/fh_icm_api:$PATH"
export PYTHONPATH="$HOME/LOMtest/STM32Workspace_1220/src/tools/python:$HOME/LOMtest/STM32Workspace_1220/src/tools/wuBase-python"
python $HOME/LOMtest/fh_server/scripts/mb_off.py -w 0
#python3 $HOME/LOMtest/fh_server/scripts/wp_off.py 
python3 $HOME/LOMtest/fh_server/scripts/wp_on.py
python $HOME/LOMtest/fh_server/scripts/mb_on.py -w 0
python3 $HOME/LOMtest/fh_server/scripts/icm_status.py

python3 $HOME/LOMtest/fh_server/scripts/icm_fpga_reboot.py -w 0 -i 2
python3 $HOME/LOMtest/fh_server/scripts/mcu_reset.py -w 0
python3 $HOME/LOMtest/fh_server/scripts/pmt_hv_enable.py -w 0
python3 $HOME/LOMtest/fh_server/scripts/icm_status.py
export LOMTEST="$HOME/LOMtest/STM32Workspace_1220/src/tools/python/LOMTest"
export TEST="/home/icehap/LOMtest/lom-test-data"
