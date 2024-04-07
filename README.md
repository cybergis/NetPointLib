# Scalable map matching and network point statistics library powered by CyberGIS

This repo contains the code and the notebooks for the on-going work about using HPC to perform map matching to network and analysis.

## Overview
Network-constrained phenomena, such as traffic accidents, crime incidents, and facility locations, are ubiquitous in urban environments. Analyzing the spatial patterns of these point events in network space is crucial for understanding the underlying processes and informing decision-making. However, the computational complexity of network-based analysis often poses significant challenges, especially when dealing with large-scale datasets or computationally intensive algorithms. To address these challenges, we present NetPointLib, a high-performance library for efficient processing and analysis of point event data in network space. NetPointLib leverages the power of the HPC virtual roger and the CyberGIS-Jupyter platform to provide a scalable and user-friendly environment for network point pattern analysis. Additionally, the library incorporates state-of-the-art implementation of multiple algorithms, such as the network local K-function and network scan statistics, enabling researchers and practitioners to uncover spatial patterns and clusters in network-constrained data.

## Requirements

### Python
- [`JPype`](https://jpype.readthedocs.io/en/latest/)
- [`numpy`](https://numpy.org/devdocs/)
- [`geopandas`](https://geopandas.org/en/stable/)
- [`pandas`](https://pandas.pydata.org/)
- [`anytree`](https://anytree.readthedocs.io/en/latest/)
- [`osmnx`](https://osmnx.readthedocs.io/en/stable/user-reference.html)

### Notebook
- [`matplotlib`](https://matplotlib.org/)


## To Start the process

```
$ python multiMapMatching_cropped.py 'place' 'data_path'  'output_path' 'x_column' 'y_column' 'number_of_processes'
```


For example:

```
$ python multiMapMatching_cropped.py Chicago Crimes_-_366.csv  Chicago_1year_nogroup Longitude Latitude 10
```

Matches the 
