# Molecular Graph Processing and Fingerprinting System

This repository contains my **Python-based molecular graph processing tool**, designed to model atomic structures as graphs and analyze chemical compounds. Built as a personal project, it features a **2048-bit fingerprinting algorithm** using **SHA-256 hashing** for substructure detection and functional group classification.


## Project Structur
molecule.py - main function, where graph stuctures and their functions are defined
molTesting.py - unit testing with PyTest
sdf.zip - data file
personalTesting folder - even more testing, less structured. Includes images of vizualized graph molecules, bonds, and atoms. 

## Overview
- **Duration**: June 2024 – August 2024  
- **Goal**: Develop a system to represent molecules as graphs, detect structural patterns, and classify compounds by functional groups.  
- **Key Results**:  
  - Substructure search via isomorphic graph matching.  
  - Unique 2048-bit fingerprints for molecular patterns.  
  - Feature extraction for compound classification.  

## Features
- **Graph Modeling**: Uses **NetworkX** to convert atomic structures into graphs for substructure analysis.  
- **Fingerprinting**: Generates **2048-bit fingerprints** with SHA-256 hashing to uniquely identify molecular patterns.  
- **Feature Extraction**: Identifies functional groups (e.g., hydroxyl, carbonyl) for classification tasks.  

## Tech Stack
- **Languages**: Python  
- **Libraries**: NetworkX, NumPy, hashlib (SHA-256), Pandas  
- **Tools**: Git, Command Line, Jupyter Notebook  
