## Master Thesis Summary
Mathematical Modeling for Short Fiber Reinforced Composites

## Goal:
Use mathematical modeling + 3D FEM to predict stiffness of short fiber
reinforced composites based on fiber orientation tensor + fiber volume fraction.

## Chapter 1:
Intro to composite materials + review of Jonathan et al. (2018).
Motivation: extend their work by including fiber volume fraction effects.

## Chapter 2:
Core mathematical tools + numerical methods (tensors, variational forms, FEM).
Clear examples + figures to explain concepts.

## Chapter 3 (Main Contribution):
Build a step-by-step computational pipeline:
- Input: fiber orientation tensor + fiber volume fraction
- Compute: stiffness tensor using 3D FEM + interpolation
Completed an open topic from Jonathan et al. by integrating volume fraction.
Includes detailed algorithm + visual workflow.

## Chapter 4:
Python implementation using Fraunhofer ITWM real data.
Compared prismatic vs. tetrahedral FEM:
- accuracy
- - computational efficiency

## Outcome:
A reproducible, mathematically grounded framework for stiffness prediction
in short fiber composites, implemented and validated on real industrial data.

## Code & Dataset:
- dataset_fraunhofer/  -> Fraunhofer ITWM dataset
- thesis_python/       -> All code used in the thesis
                          Written independently by the author in 2022, without any AI assistance.
