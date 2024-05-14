class Atom:

    def __init__(
        self,
        symbol: str,
        closures: dict[int, str] | None = None,
        charge=0,
        chirality="",
        label=0,
        nprotons=None,
    ):

        self.symbol: str = symbol
        """ List of ids of closures to other atoms"""
        self.closures: dict[int, str] = closures if closures is not None else {}
        self.charge: int = charge
        self.nprotons: int | None = nprotons
        # set to True if proton count has been set manually and need not to be overriden
        self.explicit_protons = False
        self.label: str = label
        """ SMILES chirality is @ or @@"""
        self.chirality = chirality

    def __str__(self) -> str:
        return f"{self.symbol}"

    def __repr__(self) -> str:
        t = f"'{self.symbol}'"
        if self.closures:
            t += f", {self.closures}"
        if self.charge:
            t += f", charge={self.charge}"
        if self.label:
            t += f", label={self.label}"
        return f"Atom({t})"

    def add_closure(self, closure_id: int, bond: str = "-") -> None:
        self.closures[closure_id] = bond

    def smiles_repr(self) -> str:
        t = self.symbol
        if self.closures:
            t += "".join(
                [
                    (bond if bond != "-" else "") + (f"{id}" if id < 10 else f"%{id}")
                    for id, bond in self.closures.items()
                ]
            )
        if self.nprotons:
            t += "H" + f"{self.nprotons}"
        if self.charge:
            t += f"{self.charge}"
        if self.charge or self.nprotons:
            t = f"[{t}]"
        return t

    def __hash__(self):
        return id(self.label)

    def __eq__(self, other):
        return id(self) == id(other)

    def __lt__(self, other):
        return id(self) < id(other)

    @property
    # FIXME to be removed
    def auto(self) -> bool:
        return self.nprotons is None
