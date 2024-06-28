# Option Pricing Models Repository

This repository contains MATLAB functions for pricing options using the Binomial model approach. It includes two main scripts: `Binomial.m` and `BinomialTree.m`, each serving distinct roles in the computation of option prices.

## Files Description

### 1. `Binomial.m`

The `Binomial.m` function is designed to price both European and American options (calls and puts) using various binomial tree methodologies such as the Equal Probabilities (EQP), Leisen-Reimer (LR), Cox-Ross-Rubinstein (CRR), and Tian's model. 

### 2. `BinomialTree.m`

The `BinomialTree.m` function is designed to build and evaluate a binomial tree for option pricing. It calculates the price of an option (European and American, call and put) by constructing the tree from the terminal nodes backward to the initial node, considering the possibility of early exercise for American options.
