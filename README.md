# IsingOccu
This is a working repo of Yunyi's MSc project, understanding coexistence of Muntjacs, Pheasants(SE China) and Foxes(N Wisc) using Ising Bilayer model. 
IsingOccu is a developing statistical tool to fit an Ising model with long range interaction using data with imperfect detections.
## Introduction
### Ising Model
Ising model is one of the most popular statistical mechanics models in the 20th, which models the ferromagnetism. This model use discrete variable to represent the magnetic moments (states) of atomic *spins*. Spins are usually on a graph, for example a lattice, which is a 2D Ising model. There also can exist interaction between two spins. If the interaction makes spin tending to have the same state, it is called *ferromagnetic* coupling while if not, it is called *anti-ferromagnetic* coupling. The system can be put into an *external field* which makes some spins favor certain state. It is a special case of general Probabilistic Graphical Models ([PGM](https://en.wikipedia.org/wiki/Graphical_model)) which models symmetric interactions with discrete state space.

In a 2D and 1D case, Ising model is exactly solvable. In late 20th people also studied Ising bilayer, which is the simplest setting of 3D Ising model. The local interlayer coupling can be either ferromagnetic or anti-ferromagnetic, as long as it is symmetric. Multiple approximations (e.g. Mean Field, Effective Field, Renormalizing Group, Green Function, etc.) can be made to calculate expected value of a certain physical quantity. For instance, inner energy or magnetic strength.

In statistics, the revers problem of Ising model (knowing the configuration of spins, infer the interactions) emerged in recent decades after people can manipulate the environment of lots of spin-like systems. In a simple spatial statistical setting, it is called Autologistic Model. This repo will generally follow the "single layer" autologistic model setting.

### Relation with Ecology
It is nature to consider species on the landscape with dispersion among sites to be Ising layer in external field. Dispersion or any kind of spatial relationship will introduce ferromagnetic coupling within species. If there is another species competing, another Ising layer can be introduced. Competition tend to cause partitioning, which can be view as anti-ferromagnetic interlayer coupling. Extend to three, or even four species interaction is natural as long as their interactions are symmetric (i.e. mutualism or competition) but more hard to calculate and fit a model.

The difference between a ecological setting and usual physical setting is external field. As mentioned before, 2 layer Ising model or Ising model with different types of coupling distance were intensely studied, but usually in 0 external field or homogeneous external field setting. In this kind of settings, system can be highly symmetric, e.g. Bethe approximation can be used since translational symmetry. However in a ecology setting, environment is hard to be a 0 external field (which means a totally neutral environment) or homogeneous for every species (competitional neutral). 

It is a simplified even over-simplified model of competing on landscape (e.g. the symmetric interaction), but it is not trivial, mathematically. For Ising bilayer models with slightly complex setting, MCMC is still the most powerful way to study it among physicists. Though this project is caring about different values (not critical temperature for me obviously), the behavior may be still not intuitive. MCMC and various kinds of approximations will be tried to gain insight into the model itself.

## Setting of Models
### Ising Bilayer with Inhomogeneous External Field
Just as said in the title. The competition coexistence problem is modeled as a Ising bilayer with ferromagnetic, exponential decay intralayer coupling and anti-ferromagnetic, local interlayer coupling. The whole bilayer were put into an inhomogeneous external field. One special setting here is, the external field is not only inhomogeneous through out the space, but also different between layers which make sense because different species react differently with environments. This is also how a symmetric interaction can cause asymmertric outcomes of frequency (the proportion of sites occupied).

### Reverse Problem and Statistical Model
In real animal ecology situation. One biggest problem of data is that, the detections are imperfect, i.e. data is biased to negative (absence) side because, especially foxes are shy. We have to estimate the probability of seeing the species given the species exist, using repeat sampling data. That is, we use Ising model as the "zero-inflate" distribution, and detections are Bernoulli random variables. By linking all parameters of external field and detection probabilities to environments, can we reduce the number of parameters to be fitted.

## Treatments
### Mean Field Approximation in near Neutral Environment
Mean field approximation is always the first to be tried to obtain simple property of the system. Here I use it under the near neutral setting, in which I can Taylor expand the meanfield equations near 0 and solve the linear system.

## Primary Results
### Mean Field Approximation in near Neutral Environment
![eq1](https://github.com/YunyiShen/IsingOccu/blob/master/Equations/MF.png)
In which hs are the mean environmental reactions of sigma and tau species, JNs are the cumulated intralayer coupling strength, Jsigmatau is the intralayer coupling strength, in a competition setting, it is negative.

It makes qualititive sense when the environment is nice enough, increasing dispersion will increase frequency but if its not good, increasing dispersion maybe not good (fugitive paradox). The asymmertric outcomes of frequency is because the external field though couplings are symmetric. This also has good ecological sense because niche overlap can be symmetric while the width and peak values are not, which causes some species to be stronger. This model cannot predict the frequency distribution in a total neutral environment which maybe a weak part of the theory. 


