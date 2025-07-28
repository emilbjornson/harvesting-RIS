Energy Harvesting Maximization for Reconfigurable Intelligent Surfaces Using Amplitude Measurements
==================

This is a code package related to the following scientific article:

Morteza Tavana, Meysam Masoudi, Emil Björnson, “[Energy Harvesting Maximization for Reconfigurable Intelligent Surfaces Using Amplitude Measurements](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10356096),” IEEE Transactions on Communications, vol. 72, no. 4, pp. 2201-2215, April 2024.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Energy harvesting can enable a reconfigurable intelligent surface (RIS) to self-sustain its operations without relying on external power sources. In this paper, we consider the problem of energy harvesting for RISs in the absence of coordination with the ambient RF source. We propose a series of sequential phase-alignment algorithms that maximize the received power based on only power measurements. We prove the convergence of the proposed algorithm to the optimal value for the noiseless scenario. However, for the noisy scenario, we propose a linear least squares estimator. We prove that within the class of linear estimators, the optimal set of measurement phases are equally-spaced phases. To evaluate the performance of the proposed method, we introduce a random phase update algorithm as a benchmark. Our simulation results show that the proposed algorithms outperform the random phase update method in terms of achieved power after convergence while requiring fewer measurements per phase update. Using simulations, we show that in a noiseless scenario with a discrete set of possible phase shifts for the RIS elements, the proposed method is suboptimal, achieving a higher value than the random algorithm but not exactly the maximum feasible value that we obtained by exhaustive search.


## Content of Code Package

The article contains 8 simulation figures, numbered 3a, 3b, 3c, 6, 7, 8, 9, and 10. For most figures, there is a script called simulationFigureX.m that can be run to reproduce Figure X.
An exception is Figure 10, which is reproduced by first running RIS_Channel.m to pregenerate channels, and then running simulationFigure10.m to generate the figure.

See each file for further documentation.


## Acknowledgements

The paper was supported by Digital Futures.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
