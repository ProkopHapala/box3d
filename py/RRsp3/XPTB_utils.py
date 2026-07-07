"""Utility functions for XPTB/RRsp3 molecular packing and setup.

Provides: pack_molecules_contiguous, make_h2o_geometry, masses_from_elems,
          perturb_state, load_xyz
"""

import os
import sys
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_PARENT_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..'))
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Mass lookup from elements.py if available, otherwise fallback dict
_MASS_TABLE = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'Ne': 20.18,
    'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09, 'P': 30.97,
    'S': 32.07, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.10, 'Ca': 40.08,
    'Fe': 55.85, 'Cu': 63.55, 'Zn': 65.41, 'Br': 79.90, 'I': 126.9,
    'X': 0.0,
}


def _build_mass_table():
    """Try to load full mass table from elements.py; fallback to _MASS_TABLE."""
    try:
        from elements import ELEMENTS, index_mass, index_symbol
        tbl = {}
        for row in ELEMENTS:
            sym = str(row[index_symbol])
            tbl[sym] = float(row[index_mass])
        tbl['X'] = 0.0
        return tbl
    except Exception:
        return dict(_MASS_TABLE)


_FULL_MASS = _build_mass_table()


def masses_from_elems(elems):
    """Return float32 array of atomic masses for given element symbols."""
    return np.array([_FULL_MASS.get(str(e).strip(), 0.0) for e in elems], dtype=np.float32)


def make_h2o_geometry(add_angle=False):
    """Return (pos, bonds) for a water molecule.

    Atom order: 0=O, 1=H, 2=H
    O at origin, H's at ~0.96 A with ~104.5 deg angle.

    If add_angle=True, bonds includes an angle "bond" between the two H's
    (for topology testing). Returns (pos, bonds) where bonds is a list of (i,j).
    """
    r_oh = 0.96
    angle = np.deg2rad(104.5)
    o = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    h1 = np.array([r_oh * np.cos(np.pi / 2 - angle / 2), r_oh * np.sin(np.pi / 2 - angle / 2), 0.0], dtype=np.float32)
    h2 = np.array([r_oh * np.cos(np.pi / 2 + angle / 2), r_oh * np.sin(np.pi / 2 + angle / 2), 0.0], dtype=np.float32)
    pos = np.stack([o, h1, h2], axis=0)
    bonds = [(0, 1), (0, 2)]
    if add_angle:
        bonds.append((1, 2))
    return pos, bonds


def pack_molecules_contiguous(mols, group_size=64, nodes_first=True, pad_to_group=True):
    """Pack molecules into contiguous groups of size group_size.

    Each molecule gets its own group. Within each group:
      - If nodes_first=True, node atoms (degree > 1) are placed first, then caps.
      - Remaining slots are padded with 'X' elements (if pad_to_group=True).

    Args:
        mols: list of dicts with keys 'elems', 'pos', 'bonds', 'nnode'
        group_size: atoms per group (64 for RRsp3)
        nodes_first: reorder nodes first within each group
        pad_to_group: pad each group to group_size

    Returns:
        dict with keys:
            natoms_total: int, total atom count (= n_mols * group_size if pad_to_group)
            pos: (natoms_total, 3) float32
            elems: list of element symbols (with 'X' for padding)
            is_padding: (natoms_total,) bool array
            bonds: list of (i, j) global bond indices
            nnode_group: (n_mols,) int32 array of nnode per group
    """
    group_size = int(group_size)
    nmol = len(mols)
    all_pos = []
    all_elems = []
    all_bonds = []
    all_is_pad = []
    nnode_group = np.zeros((nmol,), dtype=np.int32)

    for imol, mol in enumerate(mols):
        elems = list(mol['elems'])
        pos = np.asarray(mol['pos'], dtype=np.float32)
        bonds = list(mol['bonds'])
        n = len(elems)
        base = imol * group_size

        # Determine nodes: degree > 1
        deg = np.zeros((n,), dtype=np.int32)
        for (i, j) in bonds:
            deg[int(i)] += 1
            deg[int(j)] += 1

        if 'nnode' in mol and mol['nnode'] is not None:
            nnode = int(mol['nnode'])
        else:
            nnode = int(np.sum(deg > 1))

        nnode_group[imol] = nnode

        if nodes_first and nnode > 0 and nnode < n:
            # Reorder: nodes first, then caps
            is_node = deg > 1
            if int(np.sum(is_node)) != nnode:
                # Fallback: first nnode atoms are nodes
                is_node = np.zeros((n,), dtype=bool)
                is_node[:nnode] = True
            order = np.concatenate([np.nonzero(is_node)[0], np.nonzero(~is_node)[0]]).astype(np.int32)
            perm_inv = np.empty((n,), dtype=np.int32)
            perm_inv[order] = np.arange(n, dtype=np.int32)
            elems = [elems[i] for i in order]
            pos = pos[order, :].copy()
            bonds = [(int(perm_inv[int(i)]), int(perm_inv[int(j)])) for (i, j) in bonds]
        elif nodes_first and nnode == n:
            pass  # all atoms are nodes, no reordering needed

        # Pad to group_size
        npad = group_size - n
        if pad_to_group and npad > 0:
            pad_pos = np.zeros((npad, 3), dtype=np.float32)
            pos = np.concatenate([pos, pad_pos], axis=0)
            elems = elems + ['X'] * npad
            is_pad = np.concatenate([np.zeros((n,), dtype=bool), np.ones((npad,), dtype=bool)])
        elif pad_to_group and npad == 0:
            is_pad = np.zeros((n,), dtype=bool)
        else:
            is_pad = np.zeros((n,), dtype=bool)

        # Remap bonds to global indices
        for (i, j) in bonds:
            all_bonds.append((int(i) + base, int(j) + base))

        all_pos.append(pos)
        all_elems.extend(elems)
        all_is_pad.append(is_pad)

    pos_all = np.concatenate(all_pos, axis=0) if all_pos else np.zeros((0, 3), dtype=np.float32)
    is_pad_all = np.concatenate(all_is_pad, axis=0) if all_is_pad else np.zeros((0,), dtype=bool)

    return {
        'natoms_total': int(pos_all.shape[0]),
        'pos': pos_all,
        'elems': all_elems,
        'is_padding': is_pad_all,
        'bonds': all_bonds,
        'nnode_group': nnode_group,
    }


def perturb_state(pos, quat, pos_scale=0.05, rot_scale=0.0, rng=None):
    """Perturb atomic positions and quaternion orientations.

    Args:
        pos: (N, 3) or (N, 4) float32 positions
        quat: (N, 4) float32 quaternions [x, y, z, w]
        pos_scale: std dev of positional noise
        rot_scale: std dev of rotational noise (radians)
        rng: np.random.Generator

    Returns:
        (pos_perturbed, quat_perturbed) both float32
    """
    if rng is None:
        rng = np.random.default_rng(0)
    pos = np.asarray(pos, dtype=np.float32).copy()
    quat = np.asarray(quat, dtype=np.float32).copy()
    n = pos.shape[0]
    had_4 = (pos.shape[1] == 4)

    pos[:, :3] += rng.normal(0.0, float(pos_scale), size=(n, 3)).astype(np.float32)

    if rot_scale > 0.0:
        for i in range(n):
            axis = rng.normal(0.0, 1.0, size=3).astype(np.float32)
            norm = float(np.linalg.norm(axis))
            if norm < 1e-12:
                continue
            axis /= norm
            angle = rng.normal(0.0, float(rot_scale))
            half = angle * 0.5
            qrot = np.array([axis[0] * np.sin(half), axis[1] * np.sin(half), axis[2] * np.sin(half), np.cos(half)], dtype=np.float32)
            x1, y1, z1, w1 = qrot
            x2, y2, z2, w2 = quat[i]
            quat[i] = np.array([
                w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
                w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            ], dtype=np.float32)

    return pos, quat


def load_xyz(path):
    """Load an XYZ file.

    Returns (elems, xyz, quat) where:
        elems: list of element symbol strings
        xyz: (N, 3) float32 positions
        quat: (N, 4) float32 quaternions (all identity [0,0,0,1])
    """
    elems = []
    coords = []
    with open(path, 'r') as f:
        n = int(f.readline().strip())
        f.readline()  # comment line
        for _ in range(n):
            parts = f.readline().split()
            elems.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    xyz = np.array(coords, dtype=np.float32)
    quat = np.zeros((n, 4), dtype=np.float32)
    quat[:, 3] = 1.0
    return elems, xyz, quat
