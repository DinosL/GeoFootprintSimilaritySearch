# Similarity Search based on Geo-footprints

### General

Our implementation uses the [Superliminal RTree](https://superliminal.com/sources/#C_Code). 


### Compile

Compile using ``` g++ -O3 <program name>``` where ```<program name>``` can be either ```main_geo_footprints.cpp``` for the methods with indexing or ```normComputation.cpp``` for the methods without indexing.

### Parameters

To run the program use ```./a.out <parameters> <data_file> <query_file>```. 

Parameters for ```main_geo_footprints.cpp```:

| Parameters | README |
| ------ | ------ |
| -n | Iterative method.|
| -b | Batch search method.|
| -s | User Centric method.|
| -k | Number of similar users to be returned. Should be used together with one of the above versions (see example). |

**Example:** ./a.out -k 3 -b indexedUsers_part1.csv query_part1.csv

Parameters for ```normComputation.cpp```:

| Parameters | README |
| ------ | ------ |
| -s | Similarity (Algorithm 2 in the paper). In this method, user and query norms are pre-computed|
| -n | Similarity & Norm computation method (Algorithm 3 in the paper). Here, norms and similarity are computed simultaneously.|
| -p | Join-based similarity computation (Algorithm 4 in the paper)|
| -k | Number of similar users to be returned. |

**Example:** ./a.out -k 3 -s indexedUsers_part1.csv query_part1.csv
