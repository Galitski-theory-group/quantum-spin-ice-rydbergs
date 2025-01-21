# How to interpret the data files?

Expectation values of various operators mentioned in Table II of the paper can be calculated using the data files herein. Data from 10 Monte Carlo runs is contained in this directory. The run number of a data file is the number before ".txt" in the name of the data file. For example "nflip_3.txt" is a data file generate from run number 3.

1. To calculate the expectation value of the operator, R, divide the number in the files named "nflip_i.txt" by 4.
2. To calculate the expectation value of the operator, H_{LR}, divide the number in the files named "HLR_i.txt" by 512.
3. To calculate the expectation value of the operator, W_{LR}^2, divide the number in the files named "W2_i.txt" by 512.