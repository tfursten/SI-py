SI-py
===========================
Self-Incompatibility Simulation (Python Version)  
Version 1.0  
Release: 9/30/2013

Author
-------
Tara Furstenau  
Arizona State University  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
[Website](http://tfursten.github.io)


Description
-----------
A spatially explicit individual-based model of a diploid plant population that reproduces according to five different models of self-incompatibility.

Requirements
------------
Python 2.7

Libraries
-----------------
  * Numpy
  * Matplotlib

Usage
------------------------
```
python SI-py.py config.txt
```
Config.txt
```

1. Title
2. Mating System 0= NSI, 1= PSI, 2=GSI, 3=SSI, 4=BSI
3. Number of Pixels
4. Starting number of individuals
5. Number of generations to run
6. Data collection rate
7. Number of markers on Chromosome
8. Number of Pollen Produced (male gametes)
9. Number of Ovules Produced (female gametes)
10. Mutation Rate
11. Deleterious Mutation Rate
12. Mean Pollen Dispersal
13. Mean Seed Dispersal
```
 
