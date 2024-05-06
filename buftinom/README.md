# Smiles to IUPAC

## Convert SMILES to IUPAC names

This package contains a tools to convert SMILES to IUPAC names.


## Features

- smiles parser
- molecule representation as a graph
- algorithms on this graph to find all structures, and decomposite complex molecule into simple ones
- convert decomposed molecule to IUPAC structured representation
- Build IUPAC name of the SMILES molecule


## Usage as CLI

```bash
python -m buftinom --help

python -m buftinom -q "CCCCC" # -> "pentane"
```

### Be carefull with special bash symbols

This will NOT work:
```bash
python -m buftinom CCCCC(=O)N
```
Because `(` and `)` are special symbols in bash,

But with escaping - will

```bash
python -m buftinom "CCCCC(=O)N"
```

## Usage as library

- Check the `buftinom/__init__.py` file for public high-level methods.
- Check the `buftinom/__main__.py` file for an example of how to use the library.

### Library usage

```python
from bufintom import smiles2iupac1

smiles2iupac1("CCCCC") #  -> "pentane"
```


### Library low-level methods example

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

# to construct the structural representation of the name
iupac_naem = iupac.construct_name()

from buftinom import iupac2str

# to convert the structured name to string representetation
iupac2str(iupac_name)
```

# Limitations

- Only selected functional groups are supported
- Spyro and Bicyclo molecules are not supported
- Stereo molecules are not supported
- Aromatic rings limited support (benzene, pyridine, only)
- Esters have limites support. Some cases with multiple ester groups are not supported.

## Supported functional groups

Full list is due to be extended.
See `buftinom/fungroups.py:get_matchers` to get the full list.

Currently supported:

- ester
- acid_amide
- acid
- nitro
- oxy
- amino
- alco
- aldehyde
- ketone
- amine
- imine
- nitrile
- brom
- chlor

## Contribute new functional group support

`funcgroups.py` contains a list of functional groups and their matchers.

To add a new functional group, you need to create a new matcher function in `funcgroups.py` and add it to the `get_matchers` function.

Example of a matcher function:

```python
def nitro_matcher(mol: Molecule):
    """R-N(=O)O"""
    match = MatcherBuilder(mol, FunctionalGroup.NITRO)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="N", charge=+1).branch(
            match.atom(by="-", symbol="O", is_terminal=True, charge=-1),
            match.atom(by="-", symbol="C", is_root=True),
        ),
    )
```

Let's break it down:

- `MatcherBuilder` is a helper class to build a matcher
    It provides DSL to build a matcher for a functional group.
- `match.chain`
    - It's a helper function to chain multiple matchers.
    - It's used to match the whole functional group.
- `match.atom`
    - It's a helper function to match an atom.
    - It takes a symbol, charge, and other parameters to match an atom.
- `match.branch`
    - It's a helper function to match a branch.
    - It's used to match a part of the functional group that can have multiple variations.
    - For example, in the nitro functional group, the `N` atom can have two branches: `O` and `C`.
- `by` is used to match a bond type between atoms.
    - For example, `by="="` matches a double bond.
- `is_terminal` and `is_root` are used to mark the atom as a terminal or root atom.
    - Terminal atoms are atoms that don't have any other connections except ones that specified by matcher.
    - Root atoms are atoms that connects the functional group to the rest of the molecule.

```python
def amino_matcher(mol: Molecule):
    """- N -"""
    match = MatcherBuilder(mol, FunctionalGroup.AMINO)

    return match.chain(
        match.atom(symbol="C", is_root=True),
        match.atom(by="-", symbol="N"),
        match.atom(by="-", symbol="C", is_side_root=True),
        is_symmetric=True,
    )
```

- `is_side_root` is used when the C atom should be treated as a root atom for the side chain.
- `is_symmetric` is used when the functional group is symmetric.
    - It means that the functional group can be rotated and still be the same.
    - The matcher will try to match the functional group in all possible rotations.
```python

def acid_amide_matcher(mol: Molecule):
    """O = C - N"""
    match = MatcherBuilder(mol, FunctionalGroup.AMIDE)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="N", is_terminal=True),
        is_flex_root=True,
    )
```

- `is_flex_root` is used when the carbon atom coud be either a part of the main functional groud or not.

I.e. in the case of acid amide

* `CCCCC(=O)N` - `pentanamide` - the `C` atom is a part of the main functional group
* `c1ccccc1C(=O)N` - `phenylamide` - the `C` atom is not a part of the main functional group, but still should be matched, and do not be treated as a methyl subchain

# References

- Base naming rules was taken from http://www.adichemistry.com/organic/basics/iupac1/organic-iupac-nomenclature.html
- In case of ambiguous and unclear rules consulted The Bible - Blue Book https://iupac.org/what-we-do/books/bluebook/