Making Cell-Free Massive MIMO Competitive With MMSE Processing and Centralized Implementation
==================

This is a code package is related to the follow scientific article:

Emil Björnson and Luca Sanguinetti, “[Making Cell-Free Massive MIMO Competitive With MMSE Processing and Centralized Implementation](https://arxiv.org/abs/1903.10611),” IEEE Transactions on Wireless Communications, vol. 19, no. 1, pp. 77-90, January 2020.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Cell-free Massive MIMO is considered as a promis- ing technology for satisfying the increasing number of users and high rate expectations in beyond-5G networks. The key idea is to let many distributed access points (APs) communicate with all users in the network, possibly by using joint coherent signal processing. The aim of this paper is to provide the first comprehensive analysis of this technology under different degrees of cooperation among the APs. Particularly, the uplink spectral efficiencies of four different cell-free implementations are analyzed, with spatially correlated fading and arbitrary linear processing. It turns out that it is possible to outperform conventional Cellular Massive MIMO and small cell networks by a wide margin, but only using global or local minimum mean-square error (MMSE) combining. This is in sharp contrast to the existing literature, which advocates for maximum-ratio combining. Also, we show that a centralized implementation with optimal MMSE processing not only maximizes the SE but largely reduces the fronthaul signaling compared to the standard distributed approach. This makes it the preferred way to operate Cell-free Massive MIMO networks. Non-linear decoding is also investigated and shown to bring negligible improvements.


## Content of Code Package

The article contains 5 simulation figures, numbered 2-6. simulationFigure2a_3.m generates Figures 2(a) and 3, simulationFigure2b.m generates Figure 2(b), simulationFigure4.m generates Figure 4, simulationFigure5.m generates Figure 5, and simulationFigure6.m generates Figure 6. The package also contains 11 Matlab functions that are used by some of the scripts.

See each file for further documentation.


## Acknowledgements

E. Björnson was supported by ELLIIT and the Wallenberg AI, Autonomous Systems and Software Program (WASP). L. Sanguinetti was supported by the University of Pisa under the PRA 2018-2019 Research Project CONCEPT.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
