### Overview
This repository contains experiments of the proposed TFS-IDS approach presented in the paper 
***Incremental Discovery of Process Models Using Trace Fragments*** 
by Daniel Schuster, Niklas Föcking, Sebastiaan J. van Zelst, and Wil M. P. van der Aalst.

Corresponding author: Daniel Schuster ([Mail](mailto:daniel.schuster@fit.fraunhofer.de?subject=github-tfs-ids-approach))


### Repository Structure
* The experiments are implemented in 
`scripts/experiments_artificial.py`.


### Experiments
We provide a Dockerfile to run the conducted experiments. The following commands execute the experiments for the event log `<log_filename>` in the directory `<log_directory>` and store the results in `<results_directory>`.
```shell
docker build -t trace_fragment_supporting_ipd .
docker run -it -e LOG_FILE=<log_filename> -e AVG_TL_TO_CUT="0.2" -v <log_directory>:/usr/data -v <results_directory>:/usr/results trace_fragment_supporting_ipd
```

### Further Results
This repository contains additional results not presented in the paper.
The results/plots are located in `experimental_results/`.