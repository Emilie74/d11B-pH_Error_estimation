# d11B-pH_Error_estimation
Estimate errors for boron derived pH values

# Requirement
## Input
### A .csv file will be given
#### Column definition (fixed format)
Report error and change nothing to the input .csv file if format does not match exactly 
All column names (the first row) should be provided in input csv file
* "Sample Name"
* "Salinity"
* "SD of Salinity"
* "Temperature"
* "SD of Temperature"
* "pKb"
* "SD of pKb"
* "d11Bsw"
* "SD of d11Bsw"
* "d11Bc"
* "SD of d11Bc"
* "pHcf"
* "SD of pHcf"
Columns "pKb", "SD of pKb","pHcf", and "SD of pHcf" are empty cells to be filled
Auto fill "0.82" into empty cells in the column "SD of pHcf" before subsequent calculations.
## Output
Won't change column names (the first row)
Directly output into the .csv file by appending columns
## Formula
### pKb
pKb = -np.log10(np.exp((-8966.9-2890.53*(S**0.5)-77.942*S+1.728*(S**1.5)-0.0996*(S**2))/T+148.0248+137.1942*(S**0.5)+1.62142*S-(24.4344+25.085*(S**0.5)+0.2474*S)*np.log(T)+0.053105*(S**0.5)*T))
### pH
pKb - np.log10((d11Bsw - d11Bcarbonate) / (alpha * d11Bcarbonate - d11Bsw + 1000 * (alpha - 1)))
### Constants
* alpha = 1.0272
## Simulation steps
1. Determine pKb
   - Take sample sets of (Salinity, Temperature) by assuming normal distribution
   - Determine avg & stdev of calculated pKb for each sample set 
2. Determine pH
   - Take sample sets of (d11Bsw, d11Bc, pKb) by assuming normal distribution
   - Determine avg & stdev of calculated pH for each sample set 
