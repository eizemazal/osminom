import pytest

from osminom.structure2name.iupac import Iupac, iupac2str
from osminom.structure2name.molecule import debug_atoms
from osminom.structure2name.molecule import Molecule
from osminom.structure2name.translate import override_lang

debug_atoms(True)


def get_name(smiles):
    mol = Molecule.from_smiles(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()
    for f in iupac.alg.functional_groups.items():
        print(f)

    iupac_name = iupac.construct_name()
    print(iupac_name)
    name = iupac2str(iupac_name)
    print(name)
    return name


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C", "methane"),
        ("CC", "ethane"),
        ("C=C", "ethene"),
        ("C#C", "ethyne"),
        #
        ("CCC", "propane"),
        ("C=CC", "propene"),
        ("C#CC", "propyne"),
        #
        ("CCCC", "butane"),
        ("C=CCC", "butene"),
        ("CCC=C", "butene"),
        ("C#CCC", "butyne"),
        ("CCC#C", "butyne"),
        ("CC=CC", "but-2-ene"),
        ("CC#CC", "but-2-yne"),
        #
        ("CCCCC", "pentane"),
        ("CC=CCC#CC", "hept-2-ene-5-yne"),
        ("CCC=CC#CC", "hept-4-ene-2-yne"),
        #
        ("C=CC#C", "butene-3-yne"),
        ("CC=C=CC", "pent-2,3-diene"),
    ],
)
def test_simple_chain_name(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(C)(C)CC", "3,3-dimethylpentane"),
        ("CC(C)CC(C)C", "2,4-dimethylpentane"),
        # Note here he found that it is hexane, ill provide numeration for clearance
        #  3 21 4   56
        ("CC(CC)C(C)CC", "3,4-dimethylhexane"),
        ("CCCCC(CC)C(C)CCC", "5-ethyl-4-methylnonane"),
    ],
)
def test_decomposed_multichain(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(C)(C)CCC", "3,3-dimethylhexane"),
        ("CCC(C)CC(C)CC", "3,5-dimethylheptane"),
        ("CCCCC(C(C)C)CCCCC", "5-prop-2-yldecane"),
        ("CCC(C)C(C(C)CC)CCC", "3,5-dimethyl-4-propylheptane"),
        ("CCC(C)C(C(C)C)C=C", "4-methyl-3-prop-2-ylhexene"),
        ("CCC(C)C(C(C)CC)C=C", "4-ethenyl-3,5-dimethylheptane"),
        ("CCC(C)C(C)CCCCC(C)C", "2,7,8-trimethyldecane"),
        ("CCCC(C)C(CC)CCC", "4-ethyl-5-methyloctane"),
        ("CC(C)CCC(CC)CC", "5-ethyl-2-methylheptane"),
        ("CCC=CCC(C)C", "6-methylhept-3-ene"),
        ("CC=CCC#CC", "hept-2-ene-5-yne"),
        ("C(C)(C)C(C)(C)", "2,3-dimethylbutane"),
        ("C(C)(C)C(C)(C)CCC", "2,3,3-trimethylhexane"),
        ("CCCCC(C(C)C(C)C)(C(C)C(C)C)CCCC", "5,5-bis(3-methylbut-2-yl)nonane"),
        ("CCCC(C(C)(C)C)CC", "3-ethyl-2,2-dimethylhexane"),
        ("CCCCC(C(C)C(C)CC)C(CC)CCC", "5-butyl-6-ethyl-3,4-dimethylnonane"),
        ("CCCCCC(C(C)C(C)CC)C(CC)CCCC", "6-(3-methylpent-2-yl)-5-ethylundecane"),
    ],
)
def test_various_molecules(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC1CCCCC1", "1-methylcyclohexane"),
        ("C1CC(C)CCC1", "1-methylcyclohexane"),
        ("C1CCCCC1C", "1-methylcyclohexane"),
        ("CC1C(C)CC(C)CC1", "1,2,4-trimethylcyclohexane"),
    ],
)
def test_cycle_names(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCO", "propanol"),
        ("OCCC", "propanol"),
        ("CC(O)C", "propane-2-ol"),
        ("CC(O)C(O)C", "butane-2,3-diol"),
        ("CCC(=O)O", "propanoic acid"),
        ("CCC(O)=O", "propanoic acid"),
        ("C(O)(=O)CC(O)=O", "propane-1,3-dioic acid"),
        ("CCCN", "propanamine"),
        ("CC(C)CCCO", "4-methylpentanol"),
        ("C1CCC(O)CC1", "cyclohexanol"),
        ("C1CCCCC1C2CCCC2", "1-cyclopentylcyclohexane"),
        ("C1CCCCC1CC2CCCC2", "1-(1-cyclopentylmethyl)cyclohexane"),
        ("CCCC(C(=O)O)CCC(=O)O", "4-carboxyheptanoic acid"),
        ("CCCC(CO)CCC(=O)O", "4-hydroxymethylheptanoic acid"),
        ("C=N", "methanimine"),
        ("C#N", "methannitrile"),
        ("C(=O)N", "methanamide"),
    ],
)
def test_functional_group_naming(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1", "benzene"),
        ("C1:C:C:C:C:C:1", "benzene"),
        ("C1=CC=CC=C1", "benzene"),
        ("C1CCCCC1c1ccccc1", "1-cyclohexylbenzene"),
        ("C1CCCCCC1c1ccccc1", "1-phenylcycloheptane"),
        ("C1CCCCCC1C1=CC=CC=C1", "1-phenylcycloheptane"),
        ("c1ccccc1O", "phenol"),
        ("C1:C:C:C:C:C:1O", "phenol"),
        ("C1=CC=CC=C1O", "phenol"),
        ("C1=CC=C(O)C=C1", "phenol"),
        ("c1ccccc1CCC", "1-propylbenzene"),
    ],
)
def test_aromatic(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Base molecule, others will be derived from it
        ("CCCCCCCCCCO", "decanol"),
        ("CCCC(C(=O)(O))CCCCCCO", "7-carboxydecanol"),
    ],
)
def test_groups_as_preffixes(smiles, expected):
    """Suffix
    When carbon of the functional group is part of the parent chain

    Prefix else
    """
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        #
        ("CO", "methanol"),
        ("C(=O)O", "methanoic acid"),
        ("OCCO", "ethane-1,2-diol"),
        ("C1CC1C(CO)C", "2-cyclopropylpropanol"),
        ("C1CC1CCCC1CC1", "1-(3-cyclopropylpropyl)cyclopropane"),
        #
        ("c1cc(CCO)ccc1CC", "2-(4-ethylphenyl)ethanol"),
        ("CCCOC", "1-methoxypropane"),
        ("CC(=O)CC", "butane-2-one"),
        ("OCCC1C(CO)C1", "2-(2-hydroxymethylcyclopropyl)ethanol"),
        ("C1C(CCO)C1O", "2-(2-hydroxyethyl)cyclopropanol"),
        ("C1CC1=C", "1-methylidenecyclopropane"),
        ("CC(=O)NC", "1-methylaminoethanone"),
        ("CC(=O)OC", "methyl ethanoate"),
        ("NCCN", "ethane-1,2-diamine"),
        ("c1ccccc1N", "aminobenzene"),
        ("C1CCCCC1N", "cyclohexanamine"),
        ("C1CCCC(=O)CC(C#N)CCC1", "cyclodecanenitrile-3-one"),
        ("C1CC(C(=O)O)CCC1", "cyclohexanoic acid"),
        ("OC(=O)CC(C(=O))C(=O)O", "2-formylbutane-1,4-dioic acid"),
        ("C1CCCCC1=O", "cyclohexanone"),
        ("c1ccncc1", "pyridine"),
        ("c1ccncc(O)1", "pyrid-3-ol"),
        ("c1ccnc(N)c1", "2-aminopyridine"),
        ("c1ccccc1C(=O)N", "phenylamide"),
    ],
)
def test_more_namings(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCOCCC", "1-propoxypropane"),
        ("CCCNCCC", "1-propylaminopropane"),
        ("CCCNCCCNCCC", "1,3-dipropylaminopropane"),
        ("CCCCNCCCNCCC", "1-(3-propylaminopropylamino)butane"),
        ("OCCOCCO", "2-(2-hydroxyethoxy)ethanol"),
        ("CNC", "1-methylaminomethane"),
    ],
)
def test_splitted_by_func_group_molecules(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC=O", "propanal"),
        ("CC(=O)CC", "butane-2-one"),
    ],
)
def test_ketonealdehyde(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCC(O)O", "butane-1,1-diol"),
    ],
)
def test_many_func_groups(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCBr", "brompropane"),
        ("CC(Br)C", "2-brompropane"),
        ("C(Br)C(Br)C", "1,2-dibrompropane"),
        ("CCCCl", "chlorpropane"),
        ("CC(Cl)C", "2-chlorpropane"),
        ("C(Cl)C(Cl)C", "1,2-dichlorpropane"),
        ("CCC[N+](=O)[O-]", "nitropropane"),
    ],
)
def test_halo_groups(smiles, expected):
    assert get_name(smiles) == expected


@pytest.mark.parametrize(
    "smiles,en,ru",
    [
        ("CCCC", "butane", "бутан"),
        ("OCCCO", "propane-1,3-diol", "пропан-1,3-диол"),
    ],
)
def test_translate(smiles, en, ru):
    assert get_name(smiles) == en

    with override_lang("ru_RU"):
        assert get_name(smiles) == ru


@pytest.mark.parametrize(
    "smiles",
    [
        "C1CCCCOCCC1",
        "C1CCCCNCCC1",
    ],
)
def test_unsupported_hetero_cycles(smiles):
    with pytest.raises(NotImplementedError):
        mol = Molecule.from_smiles(smiles)
        print(iupac2str(Iupac(mol).construct_name()))
