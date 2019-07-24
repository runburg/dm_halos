# MADHAT
Model-Agnostic Dark Halo Analysis Tool

MADHAT is a numerical tool which implements a Fermi-LAT data-driven and model-independent analysis of gamma-ray emission due to dark matter annihilation in dwarf spheroidal galaxies.<br />
The analysis performed is described in https://arxiv.org/pdf/1802.03826.pdf [1]<br />
The main feature of MADHAT is that it is model agnositic in the sense it allows a user to choose J factors, thier respective errors, A list of dwarf spherodial galaxies, and DM mass and integrated energy spectrums independent of each other

## Make instructions
To make the executable run the command "make" in the MADHAT highest directory
To make the executable and remove old files run the command "make clean" in the MADHAT highest directory

In the make file the only thing that needs to be edit is line 6
```
BOOST = ../boost_1_70_0
```
Where  "boost_1_70_0" should be changed to your version number so that it can be located 

## Library Dependencies
MADHAT is dependent on the Boost Libaries and GMP Libaries. 

### Boost Libaries: [Link to website](https://www.boost.org/)

Documentation on the Boost Libaries version at the time of release can be found at this [link](https://www.boost.org/doc/libs/1_70_0/)
Installation Instructions:

Windows:
Download [boost_1_68_0.7z](https://www.boost.org/users/history/version_1_68_0.html) or [boost_1_68_0.zip](https://www.boost.org/users/history/version_1_68_0.html) and unpack it to install a complete Boost distribution

Linux and Mac OS:
Download [boost_1_68_0.tar.bz2](https://www.boost.org/users/history/version_1_68_0.html)
In the directory where you want to put the Boost installation, execute
```
tar --bzip2 -xf /path/to/boost_1_70_0.tar.bz2
```

Full installation instructions for Boost can be found at (https://www.boost.org/doc/libs/1_70_0/more/getting_started/index.html) 

### GMP Libaries: [Link to website](https://www.boost.org/)

Documentation on the Boost Libaries version at the time of release can be found at this [link](https://gmplib.org/manual/Introduction-to-GMP.html#Introduction-to-GMP)
Installation Instructions:

Full installation instructions for Boost can be found at [link](https://gmplib.org/manual/Installing-GMP.html) 

## Inputs 
The first input type is the command
```
./madhat
```
This will enter into prompts where you can choose set 0 or one of the 5 example sets used in [1]. Enter 0 to run an analysis on the custom Set_0.dat. (Input file setup described in section Configure Input Files)

1) It will then ask you if you have entered J factors corresponding to each dwarf spheroidal galaxy
2) It will then ask you if you have entered J factor errors corresponding to each J factor in 1
3) It will then ask you if for DM mass in GeV if you want to enter one
4) If you entered a DM mass then it will then ask you if for and integrated energy spectrum between 1 GeV and 100 GeV

## Outputs 
The output will be based on your inputs.

If just given a dwarf list is inputed then it will output an upper bound on the expected number of total photons arising from DM annihilation (Nbound)

and the corresponding Confidence (Beta)

If you also give it J factors then it will output Phipp in cm^3s^-1GeV^-2 (Phi)

If you also give it J factor errors then it will output Phipp erros in cm^3 s^-1 GeV^-2 (+dPhi) and (-dPhi)

Lastly if you enter the full information including DM mass in GeV and the integrated photon spectrum between 1 GeV and 100 GeV

then it will also output the annihilation cross section in cm^3 s^-1 (cs) and its errors (+dcs) and (-dcs)

Note: The cross section limit calculation assumes that the dark matter particles are their own antiparticles. If the theoretical dark matter model is Dirac instead of Majorana then divide the cross section by two. If the particle has some asymmetry between the particle and antiparticle then the factor needs to be calculated. 

## Outputs Options
There are two optional output flags built into the code that defult false 

1) on line 25 there is a bool called "fileout" if this is set to true then all input options will create and output file in the output folder.

2) on line 26 there is a bool called "printall" if this is set to true then for the first two input options it will output Nbound incrementing as an integer so that you can have Beta as a function of Nbound

## Run Options
There are 3 options to run this program

1) Set up the file Set_0.dat then use the command "./MADHAT" and it will take you through a series of prompts of what inputs are used
2) Set up a file or use Set_0.dat (The header is required), then run the command "./MADHAT Dwarf_Set Beta" i.e. "./MADHAT Set_0.dat .95" and it will only prompt about the mass and energy.
3) Set up a file or use Set_0.dat (The header is required), next set up a file of the masses and integrated energy spectrums you wish to scan over and example is mass.In, then run the command "./MADHAT Dwarf_Set Beta Input_file" i.e. "./MADHAT Set_0.dat .95 mass.in" there will be no prompts but for every line it will create a line of in a file "Output/DM_Set_Beta.out" in this example "Output/DM_Set_0_.95.out"


## Configure Input Files
To configure a run first set up the file Set_0.dat which is of the form
Note: Please only use "#" for lines that should be commented out at the top of input files
\
ID	J	dJ+	dJ-\
9	21.59	0.26	0.29\
10	21.52	0.26	0.29\
32	21.67	0.33	0.30\
36	22.25	0.37	0.62\
46	21.69	0.27	0.34\
\
Where the ID# is telling the program what dwarfs to use for the analysis and J is the Jfactor with the Log10 and dJ+ and dJ- are the corresponding errors with the log10.

The IDs can be found a the bottom of the README document

## Configure Mass Input Files
To configure a run over a set of mass and energy spectrums format the file mass.in
Note: Please only use "#" for lines that should be commented out at the top of input files
\
MASS	Energy Spectrum\
6	20.76\
10	21.51\
14	20.99\
\
Where the MASS is telling the program DM mass in GeV to use for the analysis and Energy Spectrum is the integrated photon spectrum between 1 GeV and 100 GeV

The IDs can be found a the bottom of the README document

## Dwarfs with corresponding ID number

 1--Bootes I\
 2--Bootes II\
 3--Bootes III\
 4--Canes Venatici I\
 5--Canes Venatici II\
 6--Carina\
 7--Cetus II\
 8--Columba I\
 9--Coma Berenices\
 10--Draco\
 11--Draco II\
 12--Eridanus II\
 13--Eridanus III\
 14--Fornax\
 15--Grus I\
 16--Grus II\
 17--Hercules\
 18--Horologium I\
 19--Horologium II\
 20--Hydra II\
 21--Indus II\
 22--Kim 2\
 23--Leo I\
 24--Leo II\
 25--Leo IV\
 26--Leo T\
 27--Leo V\
 28--Pegasus III\
 29--Phoenix II\
 30--Pictor I\
 31--Pisces II\
 32--Reticulum II\
 33--Reticulum III\
 34--Sagittarius II\
 35--Sculptor\
 36--Segue 1\
 37--Segue 2\
 38--Sextans\
 39--Triangulum II\
 40--Tucana II\
 41--Tucana III\
 42--Tucana IV\
 43--Tucana V\
 44--Ursa Major I\
 45--Ursa Major II\
 46--Ursa Minor\
 47--Willman I
