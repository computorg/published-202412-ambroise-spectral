# Spectral Bridges: Scalable Spectral Clustering Based on Vector Quantization

[![build and publish](https://github.com/computorg/published-202412-ambroise-spectral/actions/workflows/build.yml/badge.svg)](https://github.com/computorg/published-202412-ambroise-spectral/actions/workflows/build.yml)
[![DOI:10.57750/1gr8-bk61](https://img.shields.io/badge/DOI-10.57750/1gr8--bk61-034E79.svg)](https://doi.org/10.57750/1gr8-bk61)
[![reviews](https://img.shields.io/badge/review-report-blue)](https://github.com/computorg/published-202412-ambroise-spectral/issues/2)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/computorg/published-202412-ambroise-spectral/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/computorg/published-202412-ambroise-spectral)
[![Creative Commons License](https://i.creativecommons.org/l/by/4.0/80x15.png)](http://creativecommons.org/licenses/by/4.0/)

Authors: 

- Félix Laplante - [Université Paris-Saclay, CNRS, Univ Evry, Laboratoire de Mathématiques et Modélisation d'Evry)](http://www.math-evry.cnrs.fr/)
- [Christophe Ambroise](https://cambroise.github.io/) - [Université Paris-Saclay, CNRS, Univ Evry, Laboratoire de Mathématiques et Modélisation d'Evry)](http://www.math-evry.cnrs.fr/)

Spectral Bridges is a clustering algorithm that  builds upon the traditional k-means and spectral clustering frameworks by subdividing data into small Voronoi regions, which are subsequently merged according to a connectivity measure. Drawing inspiration from Support Vector Machine's margin concept, a non-parametric clustering approach is proposed, building an affinity margin between each pair of Voronoi regions. This approach is characterized by minimal hyperparameters and delineation of intricate, non-convex cluster structures. 

The Spectral Bridges algorithm is implemented both in 

- Python: <https://pypi.org/project/spectral-bridges>
- R: <https://github.com/cambroise/spectral-bridges-Rpackage>.

