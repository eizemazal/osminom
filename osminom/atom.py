class Atom:

    def __init__(
        self,
        symbol: str,
        conns: list[int] | None = None,
        charge=0,
        chirality="",
        label=0,
        nprotons=None,
    ):

        self.symbol: str = symbol
        """ List of ids of connections to other atoms"""
        self.conns: list[int] = conns if conns is not None else []
        self.charge: int = charge
        self.nprotons: int | None = nprotons
        self.label: str = label
        """ SMILES chirality is @ or @@"""
        self.chirality = chirality

    def __str__(self) -> str:
        return f"{self.symbol}"

    def __repr__(self) -> str:
        t = f"'{self.symbol}'"
        if self.conns:
            t += f", {self.conns}"
        if self.charge:
            t += f", charge={self.charge}"
        if self.label:
            t += f", label={self.label}"
        return f"Atom({t})"

    def add_conn(self, connid: int) -> None:
        self.conns.append(connid)

    def smiles_repr(self) -> str:
        t = self.symbol
        if self.conns:
            t += "".join([str(conn) for conn in self.conns])
        if self.nprotons:
            t += "H" + f"{self.nprotons}"
        if self.charge:
            t += f"{self.charge}"
        if self.charge or self.nprotons:
            t = f"[{t}]"
        return t
