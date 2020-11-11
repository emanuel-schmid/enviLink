# enviLink
## Access notebooks online through _binder_
**Overview:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/emanuel-schmid/enviLink/main?filepath=notebooks%2Foverview.ipynb)\
**Results:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/emanuel-schmid/enviLink/main?filepath=notebooks%2FenviLink%20results.ipynb)

Note that not all notebooks are fully operational on _binder_, some require local resources that cannot be provided online.
## Local installation and workflow execution
- download enviLink v1.0.0
- download the jar files from the release assets to bin/lib/PPS-tool-box-0.3.1-jar-with-dependencies.jar
- edit the notebooks/config.yaml file:
  - `db`: to connect to a relational database, necessary to run the _In Silico Reaction_ and _Matching_ steps 
  - `envipath`: user credentials to upload enviLink directly into envipath.org in the _Generation of Rule-Enzyme Links_ step.
- run the individual notebooks in the order following the _Overview_ notebook.
