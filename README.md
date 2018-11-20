# IsingOccu
This is a working repo of Yunyi's MSc project, understanding coexistence of Muntjacs, Pheasants(SE China) and Foxes(N Wisc) using Ising model. 
IsingOccu is a developing statistical tool to fit an Ising model with long range interaction using data with imperfect detections.

Ising model is one of the most popular statistical mechanics models in the 20th, which models the ferromagnetism. This model use discrete variable to represent the magnetic moments (states) of atomic *spins*. Spins are usually on a graph, for example a lattice, which is a 2D Ising model. There also can exist interaction between two spins. If the interaction makes spin tending to have the same state, it is called *ferromagnetic* coupling while if not, it is called *anti-ferromagnetic* coupling. The system can be put into an *outer field* which makes some spins at least favor certain state. 

People also studied Ising bilayer, which is the simplest setting of 3D Ising model. The local interlayer coupling can be either ferromagnetic or anti-ferromagnetic, as long as it is symmetric. 

It is nature to consider species on the landscape with migration between sites to be Ising layer in outer field. Migration will introduce ferromagnetic coupling. If there is another species competing, another Ising layer can be introduced. Competition tend to cause partitioning, which can be view as anti-ferromagnetic interlayer coupling. 

It is a simplified even over-simplified model of competing on landscape, but it is not trivial, mathematically. For Ising bilayer, MCMC is still the most powerful way to study it among physicists. Though I am caring about different values (not temperature for me obviously), the behavior may be still not intuitive. Thus I will use MCMC to study this system and try to obtain insight of when similar species can coexist with competition and spatial coupling.