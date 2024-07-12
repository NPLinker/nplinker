from __future__ import annotations
import json
import logging
from collections.abc import Iterator
from os import PathLike
from jsonschema import validate
from nplinker.schemas import STRAIN_MAPPINGS_SCHEMA
from .strain import Strain


logger = logging.getLogger(__name__)


class StrainCollection:
    """A collection of `Strain` objects."""

    def __init__(self) -> None:
        # the order of strains is needed for scoring part, so use a list
        self._strains: list[Strain] = []
        self._strain_dict_name: dict[str, list[Strain]] = {}

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if len(self) > 20:
            return f"StrainCollection(n={len(self)})"

        return f"StrainCollection(n={len(self)}) [" + ",".join(s.id for s in self._strains) + "]"

    def __len__(self) -> int:
        return len(self._strains)

    def __eq__(self, other) -> bool:
        if isinstance(other, StrainCollection):
            return (
                self._strains == other._strains
                and self._strain_dict_name == other._strain_dict_name
            )
        return NotImplemented

    def __add__(self, other) -> StrainCollection:
        if isinstance(other, StrainCollection):
            sc = StrainCollection()
            for strain in self._strains:
                sc.add(strain)
            for strain in other._strains:
                sc.add(strain)
            return sc
        return NotImplemented

    def __contains__(self, item: Strain) -> bool:
        """Check if the strain collection contains the given Strain object."""
        if isinstance(item, Strain):
            return item.id in self._strain_dict_name
        raise TypeError(f"Expected Strain, got {type(item)}")

    def __iter__(self) -> Iterator[Strain]:
        return iter(self._strains)

    def add(self, strain: Strain) -> None:
        """Add strain to the collection.

        If the strain already exists, merge the aliases.

        Args:
            strain: The strain to add.
        """
        if strain in self._strains:
            # only one strain object per id
            strain_ref = self._strain_dict_name[strain.id][0]
            new_aliases = [alias for alias in strain.aliases if alias not in strain_ref.aliases]
            for alias in new_aliases:
                strain_ref.add_alias(alias)
                if alias not in self._strain_dict_name:
                    self._strain_dict_name[alias] = [strain_ref]
                else:
                    self._strain_dict_name[alias].append(strain_ref)
        else:
            self._strains.append(strain)
            for name in strain.names:
                if name not in self._strain_dict_name:
                    self._strain_dict_name[name] = [strain]
                else:
                    self._strain_dict_name[name].append(strain)

    def remove(self, strain: Strain) -> None:
        """Remove a strain from the collection.

        It removes the given strain object from the collection by strain id.
        If the strain id is not found, raise `ValueError`.

        Args:
            strain: The strain to remove.

        Raises:
            ValueError: If the strain is not found in the collection.
        """
        if strain in self._strains:
            self._strains.remove(strain)
            # only one strain object per id
            strain_ref = self._strain_dict_name[strain.id][0]
            for name in strain_ref.names:
                if name in self._strain_dict_name:
                    new_strain_list = [s for s in self._strain_dict_name[name] if s.id != strain.id]
                    if not new_strain_list:
                        del self._strain_dict_name[name]
                    else:
                        self._strain_dict_name[name] = new_strain_list
        else:
            raise ValueError(f"Strain {strain} not found in the strain collection.")

    def filter(self, strain_set: set[Strain]):
        """Remove all strains that are not in `strain_set` from the strain collection.

        Args:
            strain_set: Set of strains to keep.
        """
        # note that we need to copy the list of strains, as we are modifying it
        for strain in self._strains.copy():
            if strain not in strain_set:
                self.remove(strain)

    def intersection(self, other: StrainCollection) -> StrainCollection:
        """Get the intersection of two strain collections.

        Args:
            other: The other strain collection to compare.

        Returns:
            StrainCollection object containing the strains that are in both collections.
        """
        intersection = StrainCollection()
        for strain in self:
            if strain in other:
                intersection.add(strain)
        return intersection

    def has_name(self, name: str) -> bool:
        """Check if the strain collection contains the given strain name (id or alias).

        Args:
            name: Strain name (id or alias) to check.

        Returns:
            True if the strain name is in the collection, False otherwise.
        """
        return name in self._strain_dict_name

    def lookup(self, name: str) -> list[Strain]:
        """Lookup a strain by name (id or alias).

        Args:
            name: Strain name (id or alias) to lookup.

        Returns:
            List of Strain objects with the given name.

        Raises:
            ValueError: If the strain name is not found.
        """
        if name in self._strain_dict_name:
            return self._strain_dict_name[name]
        raise ValueError(f"Strain {name} not found in the strain collection.")

    @staticmethod
    def read_json(file: str | PathLike) -> StrainCollection:
        """Read a strain mappings JSON file and return a `StrainCollection` object.

        Args:
            file: Path to the strain mappings JSON file.

        Returns:
            `StrainCollection` object.
        """
        with open(file, "r") as f:
            json_data = json.load(f)

        # validate json data
        validate(instance=json_data, schema=STRAIN_MAPPINGS_SCHEMA)

        strain_collection = StrainCollection()
        for data in json_data["strain_mappings"]:
            strain = Strain(data["strain_id"])
            for alias in data["strain_alias"]:
                strain.add_alias(alias)
            strain_collection.add(strain)
        return strain_collection

    def to_json(self, file: str | PathLike | None = None) -> str | None:
        """Convert the `StrainCollection` object to a JSON string.

        Args:
            file: Path to output JSON file. If None, return the JSON string instead.

        Returns:
            If input `file` is None, return the JSON string. Otherwise, write the JSON string to the given
            file.
        """
        data_list = [
            {"strain_id": strain.id, "strain_alias": list(strain.aliases)} for strain in self
        ]
        json_data = {"strain_mappings": data_list, "version": "1.0"}

        # validate json data
        validate(instance=json_data, schema=STRAIN_MAPPINGS_SCHEMA)

        if file is not None:
            with open(file, "w") as f:
                json.dump(json_data, f)
            return None
        return json.dumps(json_data)
