from __future__ import annotations

import numpy as np

from ..core.schema import Residue, ResidueKey

RING_SHAPE_METHOD_VERSION = "cpo-ring-shape-1"


def ring_shape_descriptors(keys: list[ResidueKey], residues: dict[ResidueKey, Residue]) -> dict[str, object]:
    """Compute rigid-transform-invariant ring shape descriptors from backbone atoms."""

    coords = []
    labels = []
    for index, key in enumerate(keys, start=1):
        for atom_name in ("N", "CA", "C"):
            atom = residues[key].atom(atom_name)
            if atom is not None:
                coords.append(atom.xyz)
                labels.append(f"{index}:{atom_name}")
    if len(coords) < 3:
        return {"method": {"name": "cyclicpepomega_ring_shape", "version": RING_SHAPE_METHOD_VERSION}, "warnings": ["fewer than three backbone atoms available"]}
    array = np.asarray(coords, dtype=float)
    centered = array - array.mean(axis=0)
    covariance = centered.T @ centered / len(centered)
    moments = sorted((float(v) for v in np.linalg.eigvalsh(covariance)), reverse=True)
    total = sum(moments) or 1.0
    normal = np.linalg.eigh(covariance)[1][:, 0]
    distances_to_plane = np.abs(centered @ normal)
    distance_matrix = []
    for i in range(len(array)):
        for j in range(i + 1, len(array)):
            distance_matrix.append(float(np.linalg.norm(array[i] - array[j])))
    rg = float(np.sqrt(np.mean(np.sum(centered ** 2, axis=1))))
    return {
        "method": {"name": "cyclicpepomega_ring_shape", "version": RING_SHAPE_METHOD_VERSION},
        "backbone_atom_count": len(coords),
        "principal_moments": moments,
        "principal_axes_variance_fraction": [moment / total for moment in moments],
        "asphericity": float(moments[0] - 0.5 * (moments[1] + moments[2])),
        "planarity_rmsd_angstrom": float(np.sqrt(np.mean(distances_to_plane ** 2))),
        "max_plane_deviation_angstrom": float(np.max(distances_to_plane)),
        "ring_compactness_radius_of_gyration_angstrom": rg,
        "distance_matrix_fingerprint": [round(value, 4) for value in sorted(distance_matrix)],
        "atom_labels": labels,
        "warnings": ["descriptor uses observed backbone atoms only"],
    }
