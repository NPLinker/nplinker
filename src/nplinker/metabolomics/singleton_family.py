from .molecular_family import MolecularFamily

class SingletonFamily(MolecularFamily):

    def __init__(self):
        super().__init__(-1)

    def __str__(self):
        return f"Singleton molecular family (id={self.id})"

