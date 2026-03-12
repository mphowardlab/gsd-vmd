"""Generate a GSD file for ethanol using JSON file from lammpsio tutorial."""

import json

import gsd.hoomd
import lammpsio
import numpy

with open("ethanol.json") as f:
    ethanol = json.load(f)

box = lammpsio.Box([-5, -5, -5], [5, 5, 5])
snap = lammpsio.Snapshot(len(ethanol["coords"]["data"]), box)

# generate particle data
unique_types = set(type_ for _, type_ in ethanol["types"]["data"])
snap.type_label = lammpsio.LabelMap(
    {idx + 1: atom_type for idx, atom_type in enumerate(unique_types)}
)
for coords, types, masses in zip(
    ethanol["coords"]["data"], ethanol["types"]["data"], ethanol["masses"]["data"]
):
    id, x, y, z = coords
    id2, type_ = types
    id3, mass = masses

    if id2 != id or id3 != id:
        raise ValueError("Data assumed to be sorted the same way")

    idx = id - 1
    snap.position[idx] = [x, y, z]
    snap.typeid[idx] = snap.type_label.inverse[type_]
    snap.mass[idx] = mass
snap.position -= numpy.sum(snap.mass[:, None] * snap.position) / numpy.sum(snap.mass)

# generate bond data
snap.bonds = lammpsio.Bonds(len(ethanol["bonds"]["data"]))
unique_bond_types = set(type_ for type_, _, _ in ethanol["bonds"]["data"])
snap.bonds.type_label = lammpsio.LabelMap(
    {idx + 1: bond_type for idx, bond_type in enumerate(unique_bond_types)}
)
for idx, bond in enumerate(ethanol["bonds"]["data"]):
    type_, i, j = bond
    snap.bonds.typeid[idx] = snap.bonds.type_label.inverse[type_]
    snap.bonds.members[idx] = [i, j]

# convert to GSD frame and write out
frame = snap.to_hoomd_gsd()
with gsd.hoomd.open("ethanol.gsd", "w") as f:
    f.append(frame)
