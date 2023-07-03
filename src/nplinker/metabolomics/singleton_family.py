import uuid
from .molecular_family import MolecularFamily


class SingletonFamily(MolecularFamily):

    def __init__(self):
        super().__init__("singleton-" + str(uuid.uuid4()))

    def __str__(self):
        return f"Singleton molecular family (id={self.id})"
