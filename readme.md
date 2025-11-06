## PyPitzer

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![PyPi license](https://badgen.net/pypi/license/pip/)](https://pypi.org/project/pip/)

[![DOI](https://zenodo.org/badge/636408216.svg)](https://zenodo.org/badge/latestdoi/636408216)

![Python Version](https://img.shields.io/badge/Python-3.12-blue.svg)

## Introduction

This package is for low-temperature themodynamic modeling with the accurate Pitzer model
and its application in quantification fluid inclsuion data and calculation of uncertainties propagated
in this modeling / quantification process.

To quantify fluid inclusion data, you the following necessary data:

* microthemometric data
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
from src.Pitzer.models import FluidPitzer
from src.Pitzer.methods import binary_query

species = {
    'Na+': 1, # always be 1 if Na is the internal standard
    'K+': 2 ,  # K/Na = 2
    # 'Mg+2':0
}

# create a fluid object with information from microthemometric and LA-ICP-MS data
fluid = FluidPitzer(
    # the initial guess
    x0=(1,1),
  
    # species defined before
    species=species,
  
    # the equilibirum of the last melting solid phase
    equilibrium='KCl(s) = K+(aq) + Cl-(aq)', 
  
    # melting temperature of the last-melting solid phase, [°C]
    t = 25,
)

result = fluid.optimize()
# result = fluid.get_b((1,1))

print(result)
```

the output will be:

```python
# Na        Cl
1.84525981, 5.53577942
```
in molality.

## Cite this package

Liu, Y., Wagner, T., & Fußwinkel, T. (2024). An integrated approach for quantifying fluid inclusion data combining microthermometry, LA-ICP-MS, and thermodynamic modeling. Chemical Geology, 644, 121863.
