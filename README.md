# pmus-benchmarking

Comparison of different literature methods of noninvasive muscle pressure estimation

### list of activities

- how to plot countour levels in heatmap (localized refining could help)
- given the low cost valley of the possible solutions (R, C, pmus) generate max and min of Pmus amplitude and start/end times (LoA).
- plot contour level just for the pixels lower than some threshold (3.0)
- how does the MIQP estimator algorithm work? does it mimic the heatmap grid search?

### list of articles, relevant concepts

- Before using the pmus_miqp.m (Marcus 2022) respiratory estimation methods, you need to set up both YALMIP and Gurobi on your machine.
- To setup YALMIP, first download YALMIP from https://yalmip.github.io/ and place its contents (YALMIP-master folder) at the root of the repository.
- Gurobi is a commercial state-of-the-art solver. If you are in academia, you are eligible to obtain a free license of Gurobi. Please visit https://www.gurobi.com/ for instructions on how to download and validate Gurobi.

### questions

- is identifiability related to the modelling or estimator (miqp, cubic, etc)? different models with distinct identifiability