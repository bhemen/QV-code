# Overview

These scripts test for equilibria in Linear Voting (LV) and Quadratic voting (QV).
See [Balancing Power in Decentralized Governance: Quadratic Voting and Information Aggregation](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4416748) for full model details.

# Model

There are $n$ voters, each voter's decision variable is the number of votes to cast, $v_i$.  We assume $v_i \in \mathbb{Z}^{+}$, i.e., $v_i$ is a non-negative integer.

We assume voter $i$ votes "correctly," i.e., with the correct sign, with probability $p_i$.  So let $s_i$ be a random variable with $\Pr [ s_i = 1] = p_i$, and $\Pr [ s_i = -1 ] = 1 -p_i$.

We say the platform reaches the "correct" decision if\
$$
    \sum_{i=1}^n s_i \cdot v_i > 0
$$
In other words, the platform reaches the correct decision if a strict majority of voters (weighted by $v_i$) votes for the "correct" outcome.

Voter $i$ receives a reward $u_i$ if the platform votes correctly and 0 otherwise.  Voter $i$'s cost is $v_i^m$ (whether or not the platform reaches the correct decision).  So voter $i$'s *expected* payoff is\ 
$$
   u_i \cdot \Pr \left[ \sum_{i=1}^n s_i \cdot v_i > 0 \right] - v_i^m
$$

We assume each voter is trying to maximize their expected payoff.

For a fixed $n$, the parameters are $p_1,\ldots,p_n$, and $u_1,\ldots,n$.  For any $n$, and any $2n$-dimensional set of parameters, we can ask **what are the equilibrium voting strategies** and **are they different under LV ($m=1$) and QV ($m = 2$)**?

# Scripts

* [QV.py](QV.py) This script has all the main functions for calculating best-responses and equilibria.  This script should not be run directly

* [search_space.py](search_space.py) This script imports [QV.py](QV.py) and searches the space of parameters for values where QV and LV have qualitatively different equilibria

* [generate_example_tex.py](generate_example_tex.py) This file takes in "interesting" parameter values and generates LaTeX tables that summarize the equilibria for these values.  This script generates the tables in the paper [The paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4416748) 
