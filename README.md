# KeNary
A Probabilistic Kernel-Based n-ary Classification Method for Sets of Observations

Code associated with the probabilistic $n$-ary classification method described in Stricklin et al. (2024).

## Introduction

In forensic science, there is oftentimes a need to classify observations into one of $n$ classes, as several pieces of evidence may be suspected to originate from a single source. In such instances, the observations should be considered in sets, rather than individually, to make an overall class determination. This code is associated with the probabilistic $n$-ary classification model, called KeNary, that allows for determining the class of a set of objects. The novelty of KeNary lies in its ability to probabilistically classify the complete set of objects at once, rather than classify each object in turn, and in that it does not require large quantities of training data. KeNary uses a kernel function to relate pairs of objects to obtain a single vector of within-class and between-class scores, and capitalizes on differences in the variability within and between these sets of scores. KeNary is inherently flexible in its ability to consider virtually any type and dimension of data, whether scalar, functional, or high-dimensional, since the kernel function at the core of the model can be modified to accommodate the data. Finally, KeNary provides a naturally probabilistic and compact multi-class alternative to current kernel-based pattern recognition methods, such as support vector machines. This method is furhter described in Stricklin et al. (2024). 

## System Requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Ventura 13.6.3, using R version 4.3.0.
Installation

To downloading and install software and packages:

    R (>= 2.14.0) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.

## Files of Interest 

See the folder **Reproduce Figures** for relevant files to recreate the figures presented in Stricklin et al. (2024). 

## Attribution and Copyright

If you use any of the KeNary framework or results in your work, please cite the following paper:

MA Stricklin, BP Weaver, JE Lee, RN Farley, RC Huber, KN Wurth, AC Aiken, KeNary Classification: A Probabilistic Kernel-Based $n$-ary Classification Method for Sets of Observations, Submitted to Technometrics.

#
Copyright 2025 for **O4858**

This program is Open-Source under the BSD-3 License.

1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

2. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

3. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

© 2025. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
