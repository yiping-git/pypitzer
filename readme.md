## PyPitzer

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![PyPi license](https://badgen.net/pypi/license/pip/)](https://pypi.org/project/pip/)

[![DOI](https://zenodo.org/badge/636408216.svg)](https://zenodo.org/badge/latestdoi/636408216)


authors: Yiping Liu, Thomas Wagner, Tobias Fußwinkel

## Introduction

This package is for low-temperature themodynamic modeling with the accurate Pitzer model
and its application in quantification fluid inclsuion data and calculation of uncertainties propagated
in this modeling / quantification process.

To quantify fluid inclusion data, you the following necessary data:

* measure microthemometry
    * accurate last solid phase (e.g. NaCl·H2O)
    * melting temperatre of this last phase (Tm)
* element/Na ratios from LA-ICP-MS analysis
    * note the elemental/Na refers to molar ratios, some conversion may be needed.

## Installation

pip install:

```python
pip install pypitzer
```

conda install:

```python
conda install pypitzer
```

## Example

After installation. Here is an example:

```python
# import the FluidPitzer class
from Pitzer.models import FluidPitzer

# aqueous species determined in LA-ICP-MS analysis
species = {
    'Na+': 1, # always be 1 if Na is the internal standard
    'K+': 2,  # K/Na = 2
}

# create a fluid object with information from microthemometric and LA-ICP-MS data
fluid = FluidPitzer(
    # the initial guess
    x0=(3, 3),
  
    # species defined before
    species=species,
  
    # the last melting solid
    solids=['KCl'],
  
    # melting temperature of the last solid, °C
    t = 25,
)

result = fluid.optimize()

print(result)
```

the output will be:

```python
# Na        Cl
1.84525981, 5.53577942
```
in molality.
## Cite this package
...