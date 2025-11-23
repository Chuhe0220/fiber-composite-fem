# ===============================================================
# fiber-composite-fem
# ===============================================================
# A mathematical and 3D FEM-based framework for predicting the
# stiffness of short fiber reinforced composites using fiber
# orientation tensors and fiber volume fractions.
#
# This repository contains the core implementation from my
# master's thesis.
# ===============================================================


# -----------------------------
# Overview
# -----------------------------
# Short fiber reinforced composites are anisotropic materials.
# Their mechanical behavior depends mainly on:
#   - fiber orientation tensor
#   - fiber volume fraction
#
# This project builds a 3D finite element + interpolation
# pipeline that computes stiffness tensors from these inputs.
#
# Highlights:
#   * Extends Jonathan et al. (2018) by incorporating
#     fiber volume fraction.
#   * Uses real industrial data from Fraunhofer ITWM.
#   * Compares prismatic vs. tetrahedral FEM in accuracy
#     and computational efficiency.
# -----------------------------


# -----------------------------
# Features
# -----------------------------
# - Compute stiffness tensors from fiber orientation
# - FEM-based interpolation across fiber states
# - Support for prismatic and tetrahedral meshes
# - Fully reproducible Python implementation
# - Includes notebooks, figures, and results
# -----------------------------


# -----------------------------
# Installation
# -----------------------------
# Clone the repository:
#   git clone https://github.com/yourname/fiber-composite-fem.git
#   cd fiber-composite-fem
#
# Install dependencies:
#   pip install -r requirements.txt
# -----------------------------


# -----------------------------
# Usage
# -----------------------------
# Run the full computation pipeline:
#   python src/run_pipeline.py
#
# Or explore the Jupyter examples:
#   notebooks/
# -----------------------------


# -----------------------------
# Project Structure
# -----------------------------
# fiber-composite-fem/
# ├── src/               # FEM + interpolation code
# ├── data/              # Fiber orientation & material data
# ├── notebooks/         # Computational experiments
# ├── figures/           # Plots and illustrations
# ├── results/           # Output stiffness tensors
# ├── requirements.txt
# └── thesis.pdf
# -----------------------------


# -----------------------------
# Reference
# -----------------------------
# Jonathan, K. S., et al.
# "Fiber orientation interpolation for the multiscale analysis
#  of short fiber reinforced composite parts" (2018).
# -----------------------------


# -----------------------------
# License
# -----------------------------
# MIT License (recommended)
# -----------------------------
