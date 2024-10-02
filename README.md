# HIWT-GSC: Iterative-Weighted Thresholding Method for Group-Sparsity-Constrained Optimization

[![License](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![MATLAB](https://img.shields.io/badge/MATLAB-blue.svg)](https://www.mathworks.com/products/matlab.html)

This repository provides the official MATLAB implementation of the **Homotopy Iterative-Weighted Thresholding** algorithm for solving group-sparsity-constrained optimization problems , as described in our paper:

**Iterative-Weighted Thresholding Method for Group-Sparsity-Constrained Optimization with Applications**  
_Author_: Lanfan Jiang, Zilin Huang, Yu Chen and Wenxing Zhu  
_Journal_: IEEE Transactions on Neural Networks and Learning Systems  
_DOI_: https://doi.org/10.1109/TNNLS.2024.3454070

> Abstract:
Taking advantage of the natural grouping structure inside data, group sparse optimization can effectively improve the efficiency and stability of high-dimensional data analysis, and it has wide applications in a variety of fields such as machine learning, signal processing, and bioinformatics. Although there has been a lot of progress, it is still a challenge to construct a group sparse-inducing function with good properties and to identify significant groups. This article aims to address the group-sparsity-constrained minimization problem. We convert the problem to an equivalent weighted $\ell_{p,q}$-norm $(p>0,0<q\leq1)$ constrained optimization model, instead of its relaxation or approximation problem. Then, by applying the proximal gradient method, a solution method with theoretical convergence analysis is developed. Moreover, based on the properties proved in the Lagrangian dual framework, the homotopy technique is employed to cope with the parameter tuning task and to ensure that the output of the proposed homotopy algorithm is an $L$-stationary point of the original problem. The proposed weighted framework, with the central idea of identifying important groups, is compatible with a wide range of support set identification strategies, which can better meet the needs of different applications and improve the robustness of the model in practice. Both simulated and real data experiments demonstrate the superiority of the proposed method in terms of group feature selection accuracy and computational efficiency. Extensive experimental results in application areas such as compressed sensing, image recognition, and classifier design show that our method has great potential in a wide range of applications.

> Keywords:
Group sparse, sparse optimization, iterative weighted thresholding, homotopy, proximal gradient, non-convex optimization.

## Table of Contents

- [Introduction](#introduction)
- [Paper and Appendix](#paper-and-appendix)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
  - [Demo Scripts](#demo-scripts)
  - [Running the Solver](#running-the-solver)
- [Citation](#citation)
- [License](#license)

## Introduction

The **HIWT-GSC** method is designed for optimizing group-sparsity-constrained problems, which arise in a variety of signal processing and machine learning applications. 
<!--
The **HIWT-GSC** repository provides a MATLAB-based implementation of the **Iterative-Weighted Thresholding (IWT)** algorithm, designed to solve group-sparsity-constrained optimization problems. This method is particularly useful in applications involving sparse recovery, compressed sensing, and signal processing.
-->
This repository includes:

- MATLAB implementation of the HIWT-GSC algorithm, a solver for group-sparsity-constrained optimization problems.
- Some demos showing how to use the HIWT_GSC for a group-sparsity-constrained problem.
- The peer-reviewed accepted version of the paper and Supplementary Material (appendix). The appendix contains all proofs of theorems, lemmas, additional discussions and experiments.

## Paper and Appendix

 Due to space limitations, the proofs of the lemmas and theorems are not presented in the main body of the paper, but are provided in the Supplementary Material, as appendix. 

You can find the peer-reviewed accepted version of the paper along with the detailed appendix in the following files:

- [Full Paper (PDF)](./paper/tnnls-hiwtgsc-fullpaper-20241002.pdf): This document contains the complete text of the paper, including both the main body and the appendix. The appendix includes all proofs of theorems and lemmas presented in the paper, as well as additional discussions and experimental results. 
The content in the appendix has also undergone peer review.
Cross-referencing are enabled for seamless navigation and enhanced readability.

We encourage readers to refer to both documents to gain a comprehensive understanding of the methodology and results discussed in our work.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/YourUsername/HIWT-GSC.git
   cd HIWT-GSC
   ```

2. Ensure you have MATLAB installed. The code has been tested on MATLAB R2016b, but it should also be compatible with any version higher than R2011a.

## Repository Structure

The repository is structured as follows:
```bash
HIWT-GSC/
│
├── code/                          # MATLAB source code
│   ├── HIWT_GSC.m                 # The outer homotopy algorithm.  
│   ├── IWT_GSC.m                  # The inner iterative algorithm, which is called by the HIWT_GSC
│   ├── demo_LeastSquares.m        # Least squares demo
│   ├── demo_LogisticRegression.m  # Logistic regression demo
│   ├── utils/                     # Utility functions
│
├── paper/
│   ├── tnnls-hiwtgsc-fullpaper-[date].pdf     # Peer-reviewed accepted version of the paper and appendex
│
├── README.md                      # Project documentation
└── LICENSE                        # License information
```

## Usage

### Demo Scripts

We provide several demo scripts to illustrate how to use the HIWT_GSC solver:

- `demo_LeastSquares.m`: Shows how to use the HIWT_GSC for a group-sparsity-constrained least squares regression problem.
- `demo_LogisticRegression.m`: Shows how to use the HIWT_GSC for a group-sparsity-constrained logistic regression problem.

To run a demo, simply execute the corresponding script in MATLAB:
```matlab
>> demo_LeastSquares
>> demo_LogisticRegression
```

### Running the Solver

The core solver is implemented in the file `HIWT_GSC.m`. It leverages a homotopy approach to gradually adjust the regularization parameter, ensuring better convergence properties when solving group-sparse-constrained optimization problems. 

You can solve group-sparsity-constrained optimization problems by calling the function with the appropriate input arguments. 

Example usage:
```matlab
result = HIWT_GSC(fun_obj,A,b,x0,options);
```

- `fun_obj`: Objective function handle 
- `A`: Sample matrix  
- `b`: Sampled data 
- `x0`: Initial solution, defaults to 0
- `options`: A structure variable determining various options of the algorithm

Refer to the comments in the function for a detailed explanation of input formats and options.


## Citation

If you use this code in your research, please cite the following paper:

```bibtex
@article{jiang-hiwt-gsc-2024,
  author = {Jiang, Lanfan and Huang, Zilin and Chen, Yu and Zhu, Wenxing},
  title = {Iterative-Weighted Thresholding Method for Group-Sparsity-Constrained Optimization with Applications},
  journal = {IEEE Transactions on Neural Networks and Learning Systems},
  year = {2024},
  note = {Early Access},
  doi = {10.1109/TNNLS.2024.3454070},
  publisher = {IEEE}
}
```

## Acknowledgments

We would like to express our sincere gratitude to the esteemed editors and the dedicated anonymous reviewers for their insightful suggestions and constructive feedback. Their contributions have been instrumental in enhancing the rigor and clarity of this paper.

We would like to express our gratitude to the authors of the paper ["Group Sparse Recovery via the $\ell_0(\ell_2)$ Penalty: Theory and Algorithm"](https://jszy.whu.edu.cn/jiaoyuling/en/lwcg/1349484/content/54896.htm) for providing the associated code. This code, which can be found [here](http://www0.cs.ucl.ac.uk/staff/b.jin/software/gpdasc.zip), has greatly supported our research on the Group-Sparsity-Constrained optimization algorithms.


## License

This project is licensed under the CC BY-NC 4.0 License - see the [LICENSE](LICENSE) file for details.
