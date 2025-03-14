# Overview

These scripts test for equilibria in Linear Voting (LV) and Quadratic voting (QV).
See [Balancing Power in Decentralized Governance: Quadratic Voting and Information Aggregation](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4416748) for full model details.

# Scripts

* [QV.py](QV.py) This script has all the main functions for calculating best-responses and equilibria.  This script should not be run directly

* [search_space.py](search_space.py) This script imports [QV.py](QV.py) and searches the space of parameters for values where QV and LV have qualitatively different equilibria

* [generate_example_tex.py](generate_example_tex.py) This file takes in "interesting" parameter values and generates LaTeX tables that summarize the equilibria for these values.  This script generates the tables in the paper [The paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4416748) 