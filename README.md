# SteinerAnts

This repository is for analyzing the design principles guiding the structure of turtle ant trail networks. This code for the paper "Better tired than lost: turtle ant trail networks favor coherence over shortest edges" in Plos Computational Biology (DOI ---).

The code for this project is  written in Python 2.7. The exact list of packages is listed in the file environment.yml. We recommend using the Anaconda python distribution, and using environment.yml to create an enviromnent that mirrors the one used for this project as follows:

    conda env create -f environment.yml

First run the following command:

    python objective_percentiles.py --write -r 10000
    
This will run 10000 random simulations using the default parameters that were used in the paper. This will create two files: objective_values.csv and objective_percentiles.csv. The latter will be used to recreate the visualizations used in the paper and perform statistical tests.

To recreate the the panels in figures 4 and 5, run the following:

    python objective_percentiles.py --plot --swarm
    
To reproduce the results of performing Friedman's test, run the followng:
    
    python objective_percentiles.py --friedman
    
To reproduce the results of performing Conover's test, run the following:

    Rscript objectives_anova.R
    
Note that this R script requires that your R installation includes the packages rcompanion, fsa, scales, stats, ARTool, PMCR, and conover.test.
    
To reproduce the results of comparing observed and optimized networks, run the following:

    python objective_percentiles.py --compare_heuristic
    
To reproduce the results of analyzing loops in random networks, first run the following:

    python multipaths.py -w -r 100000
    
This will sample from all of the loops created in observed trails. This will create a file called multipaths.csv. Then run the following command:

    python objective_percentiles.py --loops
    
To calculate correlations between pairs of metrics, run the following command:

    python network_stats.py --correlations
    
To reproduce the stats in table 1, run the following command:

    python network_stats.py --table1
