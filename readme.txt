Meltpond detection from MODIS MOD09 data:

1. Download data
2. shell skript: execute_reproj_via_shell.sh to do reprojection of hdf tiles....
3. Use read_MOD09XX_2grd_snow.py to create netcdfs
4. Use mosaicMOD09.py to create Mosaics
5. Use neural_net_applyMOD09.py  to apply ANN
6. Use mp_regrid.... to grid it to 12.5 km.....

for validation
all programms starting with validation.....

for time series analysis:
most recent: zeitreihe_5masked.py

