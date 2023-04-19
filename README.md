# transcript_feature_counter

This code can be used for transcript feature counting associated with Nonsense-mediated decay. This code uses the results of the NMDNAR program, resulting in values, with values obtained from calculations using information from this repository. 

It can be used for counting such variables:

![covariates_all_v3](https://user-images.githubusercontent.com/111967607/233170420-b02d2c7d-510a-4a88-bf3a-f323e5a0e4ed.png)

The choice of metrics was made on the basis of already available experimental data, with the addition of new metrics, which influence on NMD is not confirmed. So for example metrics (1),(3),(4),(10),(11) were taken on the basis of data from the [article](https://www.cell.com/molecular-cell/pdf/S1097-2765(19)30361-2.pdf). Other metrics were found by combinatorial enumeration of transcript characteristics that can numerically describe each transcript containing stop codon. In this way, 14 metrics were derived, and their values will be used to find significant changes in delta-PSI depending on the metric value. 

This repo can be used in two ways:

1. Counting all values and getting all consistent results
```
make all
```
2. Only get pictures based on already calculated features
```
make pic
```
