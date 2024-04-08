# Smiles to IUPAC

## Convert SMILES to IUPAC names

This package contains a command line tool to convert SMILES to IUPAC names.

```bash
python -m buftinom --help
```

## Usage as library

- Check the `buftinom/__inut__.py` file for public high-level methods.

- Check the `buftinom/__main__.py` file for an example of how to use the library.

- To use low-level check

```python
from buftinom import SmilesParser, Iupac

mols = SmilesParser().parse(smiles)

names = []

for mol in mols:
    iupac = Iupac(mol)

```

from the `Iupac` object you'll have acces to the following attributes:

```python
iupac.mol  # the molecule object
iupac.decompostion  # decomposition of the molecule
iupac.alg  # the algorithms object used to create decomposition

iupac_naem =iupac.construct_name()  # will construct the structural representation of the name

from buftinom import iupac2str

iupac2str(iupac_name)  # will convert the structured name to string representetation
```