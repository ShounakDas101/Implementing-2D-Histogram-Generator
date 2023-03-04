# Implementing-2D-Histogram-Generator
Title	:	 MadAnalysis5 - Implementing 2D Histogram Generator

Author	: 	Shounak Das

Date	:	04-03-2023

Program Code :

    Histo.cpp

Compilation:

	$ g++ Histo.cpp -o main.o

Execution:

	$./main.o input.dat distribution_index num_bins X_min X_max
  
	distribution_index (integer value)
	-----------------------
	 0   for     px
	 1   for     py
	 2   for     pz
	 3   for     E
	 4   for     pT
	 5   for     M         \
	-----------------------
	num_bins   ( integer value)
	X_min       (double value)
	X_max       (double value)
	
	example 
  
  $./main.o input.dat 1 100 -2000.87 10000.9234
