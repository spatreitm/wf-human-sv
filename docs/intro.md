## Introduction

The pipeline performs the following steps:
* Maps reads using [lra](https://github.com/ChaissonLab/LRA)
* Calls variants using [cuteSV](https://github.com/tjiangHIT/cuteSV)
* Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)
* Optionally evaluates yours calls against a truthset using [truvari](https://github.com/spiralgenetics/truvari)
