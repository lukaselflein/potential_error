# Potential Error
## Motivation
This side project aims to resolve the problem of quality of point charges.
So far, we only look at proxy metrics to measure quality: the sum of differences between constrained and unconstrained charges, inspecting those differences visually, the dispersion in the distribution of charges, and the match of human intuitions with charges.

In the original paper, the authors compared the potential calculated with their optimized point charges to the potential obtained by DFT calculations. So far, we were not able to reproduce this, as they did not publish the code they used. This project replicates this comparison.

## Usage
* mkdir data
* python crawl_calc.py
* python plot_crawl.py

Sadly, time constraints do not allow for a full-blown integration into the current workflow. At this point, the project assumes that data files are located in a folderstructure on the same hierarcical level as the project root folder.
Calling `crawl_calc.py` will start a search through `../{}_charge/1_charge_cycle` with {} ranging from 0 to 2.
All `esp.cube` files found will be compared to their respective point-charge file, assumed to be in `4_horton_cost_fuction/lnrho/` and being named `charges_-9_0.8.csv`, with -9 and 0.8 being the lnrhoref/sigma parameters of the parameter sweep.

To speed up the calculation process, the `crawl_calc.py` file is parallelized: calling it multiple times will result in performance roughly linear in the number of calls:
```
for i in `seq 1 1 20`
  do
  python crawl_calc.py > /dev/null &
  done
```
