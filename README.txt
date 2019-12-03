env: VS2015 + opencv3.1 + matlab2013
PFDTV script:
	(1): change file paths in the cpp
	(2): the script will run PAS.m to calculate PAS map during each iteration
	(3): the script will read the PCMAP.txt to despeckle
PAS.m script: 
	(1): calculate the phase asymmetry map after each iteration;
	(2): save the results in the PCMAP.txt
results：There are 200 original ultrasound images and corresponding despeckling results in the results directory. All the ultrasound images are downloaded from the public dataset.