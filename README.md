# Spectral Bridges
Félix Laplante, Christophe Ambroise
2024-12-13

[![build and
publish](https://github.com/computorg/published-202412-ambroise-spectral/actions/workflows/build.yml/badge.svg)](https://github.com/computorg/published-202412-ambroise-spectral/actions/workflows/build.yml)
[![reviews](https://img.shields.io/badge/review-report-blue.png)](https://github.com/computorg/published-202412-ambroise-spectral/issues?q=is%3Aopen+is%3Aissue+label%3Areview)
[![DOI:10.57750/1gr8-bk61](https://img.shields.io/badge/DOI-10.57750%2F1gr8--bk61-034E79.svg)](https://doi.org/10.57750/1gr8-bk61)
[![Creative Commons
License](https://i.creativecommons.org/l/by/4.0/80x15.png)](http://creativecommons.org/licenses/by/4.0/)

### Authors

- Félix Laplante (Université Paris-Saclay, CNRS, Univ Evry,)
- [Christophe Ambroise](https://cambroise.github.io/) (Université Paris-Saclay, CNRS, Univ Evry,)

### Abstract

In this paper, Spectral Bridges, a novel clustering algorithm, is
introduced. This algorithm builds upon the traditional k-means and
spectral clustering frameworks by subdividing data into small Voronoï
regions, which are subsequently merged according to a connectivity
measure. Drawing inspiration from Support Vector Machine’s margin
concept, a non-parametric clustering approach is proposed, building an
affinity margin between each pair of Voronoï regions. This approach
delineates intricate, non-convex cluster structures and is robust to
hyperparameter choice. The numerical experiments underscore Spectral
Bridges as a fast, robust, and versatile tool for clustering tasks
spanning diverse domains. Its efficacy extends to large-scale scenarios
encompassing both real-world and synthetic datasets. The Spectral Bridge
algorithm is implemented both in Python
(<https://pypi.org/project/spectral-bridges>) and R
<https://github.com/cambroise/spectral-bridges-Rpackage>).
