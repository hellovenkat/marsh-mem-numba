##Libraries needed to be installed:
+ numba
+ gdal
+ rasterio
+ numpy

##Files description:
MEM_testing_170821_DEM5.py is the original file which is single threaded.
MEM_testing_numba_DEM5.py is the file which uses multi-threading concept with Numba library to run the program.

##To run the original single-threaded program: 
From the terminal window and change to the local directory where the downloaded repository is and cd/marsh-mem-Numba.
Now just run the command, 
```
#!command

python MEM_testing_170821_DEM5.py
```

The default location for the output files is 
```
#!command

current_directory/DEM5/single_threaded
```

##To run the multi-threaded program(uses numba library): 
From the terminal window and change to the local directory where the downloaded repository is and cd/marsh-mem-Numba.

Now just run the submit.sh file which has the command,
```
#!command

python MEM_testing_numba_DEM5.py
```
 
Also the number of cores and nodes are specified in the .sh file.

The default location for the output files is 
```
#!command

current_directory/DEM5/multi_threaded
```

Note that the 'Beaufort_DEM2013_meters.tif' file and the 'MEM_testing_numba_DEM5.py' python file must be in the same folder. 
The data needed to run is extracted from 'Beaufort_DEM2013_meters.tif' by the program.
The Output files are '.tif' files, one file will be generated for each year. The default number of years is 10.
You can specify the number of years in the python file by just re-intialising the variable 'T' to the desired number of years.


