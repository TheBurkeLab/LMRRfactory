## Installation
```bash
pip install LMRRfactory
```

[![PyPI version](https://badge.fury.io/py/LMRRfactory.svg)](https://pypi.org/project/LMRRfactory/)

## How to Use
```bash
import LMRRfactory
```
To apply LMR-R to all eligible pressure-dependent reactions in your input mech for which ab initio, reaction-specific, temperature-dependent third-body efficiencies are available in the LMRRfactory internal database:
```bash
LMRRfactory.makeYAML(mechInput="your_base_model.yaml", outputPath="folder/where/you/want/the/output/file/stored")
```
To apply ab initio LMR-R as described above AND apply "generic", temperature-independent third-body efficiencies to the remaining pressure-dependent reactions which lack an ab initio set in the database:
```bash
LMRRfactory.makeYAML(mechInput="your_base_model.yaml", outputPath="folder/where/you/want/the/output/file/stored", allPdep=True)
```
To apply LMR-R only to a single reaction (e.g. H + OH (+M) <=> H2O (+M)) and leave the rest of the input mechanism unchanged:
```bash
LMRRfactory.makeYAML(mechInput="your_base_model.yaml", outputPath="folder/where/you/want/the/output/file/stored", reaction=' H + OH (+M) <=> H2O (+M)')
```
To apply LMR-R only to a single reaction AND only provide that reaction with a single third-body efficiency of your choice, e.g. HE:
```bash
LMRRfactory.makeYAML(mechInput="your_base_model.yaml", outputPath="folder/where/you/want/the/output/file/stored", reaction=' H + OH (+M) <=> H2O (+M)', collider='HE')
```

## Citation

If you use this tool, please cite:

```bibtex
@misc{LMRRfactory,
  author = {Singal,  Patrick J. and Burke,  Michael P.},
  title = {LMRRfactory: A Software Toolkit and Database for Implementing Mixture Rules at Scale},
  year = {2026},
  howpublished = {https://github.com/TheBurkeLab/LMRRfactory},
  doi = {10.5281/ZENODO.14649194},
}
```

## Acknowledgement

The authors gratefully acknowledge support from the Department of Energy Gas Phase Chemical Physics program (DE-SC0019487).
