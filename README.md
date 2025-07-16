# pmus-benchmarking

Comparison of different literature methods of noninvasive muscle pressure estimation

### list of activities

- loads ASL data in repo
- adds data processing and plotting tools
- adds marks/cycles identification tools
- try to reproduce articles behaviour in a set of 10 cycles

## list of articles, relevant concepts

- Before using the pmus_miqp.m (Marcus 2022) respiratory estimation methods, you need to set up both YALMIP and Gurobi on your machine.
- To setup YALMIP, first download YALMIP from https://yalmip.github.io/ and place its contents (YALMIP-master folder) at the root of the repository.
- Gurobi is a commercial state-of-the-art solver. If you are in academia, you are eligible to obtain a free license of Gurobi. Please visit https://www.gurobi.com/ for instructions on how to download and validate Gurobi.