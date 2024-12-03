# KeNary
A Probabilistic Kernel-Based n-ary Classification Method for Sets of Observations

Code associated with the probabilistic $n$-ary classification method described in Stricklin et al. (2024).

## Introduction

In forensic science, there is oftentimes a need to classify observations into one of $n$ classes, as several pieces of evidence may be suspected to originate from a single source. In such instances, the observations should be considered in sets, rather than individually, to make an overall class determination. This code is associated with the probabilistic $n$-ary classification model, called KeNary, that allows for determining the class of a set of objects. This method is furhter described in Stricklin et a. (2024). The novelty of KeNary lies in its ability to probabilistically classify the complete set of objects at once, rather than classify each object in turn, and in that it does not require large quantities of training data. KeNary uses a kernel function to relate pairs of objects to obtain a single vector of within-class and between-class scores, and capitalizes on differences in the variability within and between these sets of scores. KeNary is inherently flexible in its ability to consider virtually any type and dimension of data, whether scalar, functional, or high-dimensional, since the kernel function at the core of the model can be modified to accommodate the data. Finally, KeNary provides a naturally probabilistic and compact multi-class alternative to current kernel-based pattern recognition methods, such as support vector machines.

## System Requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Ventura 13.6.3, using R version 4.3.0.
Installation

To downloading and install software and packages:

    R (>= 2.14.0) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.

## Files of Interest 

Technometrics_Submit contains the various codes and files for reproducing the figures contained in Stricklin et al. (2024):

1. 1_PaintSims_Study.R will recreate the Paint study in Stricklin et al. (2024).
2. 2_SubmitFigs.R will recreate the figures in Stricklin et al. (2024).

Note that reproducing some of the figures in Stricklin et al. (2024) requires loading large .rda files (e.g, FTIRSimsResults_2sources_3test_24Sep24.rda). These files can be obtained by contacting the author at mstricklin@lanl.gov. 

## Attribution and Copyright

If you use any of the KeNary framework or results in your work, please cite the following paper:

MA Stricklin, BP Weaver, JE Lee, RN Farley, RC Huber, KN Wurth, AC Aiken, KeNary Classification: A Probabilistic Kernel-Based $n$-ary Classification Method for Sets of Observations, Submitted to Technometrics.
