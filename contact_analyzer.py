"""
contact_analyzer.py — Extract lectin-sugar contacts from PDB files using BioPython.

Replaces PLIP when it won't install. Counts:
  - H-bonds (donor-acceptor distance < 3.5 Å between sugar O/N and protein O/N)
  - Aromatic contacts (sugar C atoms within 4.5 Å of Trp/Tyr/Phe ring centroids)
  - Water bridges (water O within 3.5 Å of both sugar and protein atoms)
  - Metal coordination (Ca/Mn within 2.8 Å of sugar O)

These are APPROXIMATE counts. Manual verification in PyMOL is still required,
especially for CH-π contacts which require visual inspection of geometry.
"""

import csv
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    from Bio.PDB import PDBParser, NeighborSearch
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("WARNING: BioPython not installed. Run: pip install biopython")


# Ligand component IDs for each PDB
LIGAND_COMP_IDS = {
    "5cna": ["MMA"],         # methyl alpha-D-mannopyranoside
    "1cvn": ["MAN", "BMA"],  # mannose, beta-mannose
    "2uvo": ["NAG"],         # N-acetylglucosamine
    "2pel": ["GAL", "GLA"],  # galactose
    "1a3k": ["GAL", "NAG"],  # Galectin-3 CRD + LacNAc (Gal-β1,4-GlcNAc)
    "1lzb": ["NAG"],         # N-acetylglucosamine (chitobiose)
}

AROMATIC_RESIDUES = {"TRP", "TYR", "PHE", "HIS"}
AROMATIC_ATOMS = {
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
# Six-membered ring atoms for ring normal calculation
AROMATIC_6RING = {
    "TRP": ["CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"],  # 6-ring of indole
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],  # 5-ring
}
METAL_ELEMENTS = {"CA", "MN", "ZN", "MG"}
HB_CUTOFF = 3.5    # Å, generous cutoff for H-bond (donor-acceptor)
CH_PI_DIST_MIN = 3.0  # Å, minimum C to centroid for CH-π
CH_PI_DIST_MAX = 4.5  # Å, maximum C to centroid for CH-π
CH_PI_ANGLE_MAX = 60.0  # degrees from ring normal — within cone = CH-π
WATER_CUTOFF = 3.5  # Å, water bridge
METAL_CUTOFF = 2.8  # Å, metal coordination


def _distance(a1: Atom, a2: Atom) -> float:
    return np.linalg.norm(a1.get_vector().get_array() - a2.get_vector().get_array())


def _is_polar(atom: Atom) -> bool:
    """Is this atom a potential H-bond donor/acceptor (O or N)?"""
    elem = atom.element.strip().upper()
    return elem in ("O", "N", "S")


def _ring_centroid(residue: Residue, ring_atom_names: List[str]) -> Optional[np.ndarray]:
    """Compute centroid of aromatic ring atoms."""
    coords = []
    for name in ring_atom_names:
        if name in residue:
            coords.append(residue[name].get_vector().get_array())
    if len(coords) >= 4:
        return np.mean(coords, axis=0)
    return None


def _ring_normal(residue: Residue, ring_atom_names: List[str]) -> Optional[np.ndarray]:
    """Compute normal vector of aromatic ring plane."""
    coords = []
    for name in ring_atom_names:
        if name in residue:
            coords.append(residue[name].get_vector().get_array())
    if len(coords) >= 3:
        # Use first 3 non-collinear atoms to define plane
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        normal = np.cross(v1, v2)
        norm = np.linalg.norm(normal)
        if norm > 1e-6:
            return normal / norm
    return None


def _ch_pi_angle(carbon_pos: np.ndarray, centroid: np.ndarray,
                 ring_normal: np.ndarray) -> float:
    """Angle (degrees) between C→centroid vector and ring normal.
    0° = directly above ring (ideal CH-π). 90° = in plane (not CH-π)."""
    vec = carbon_pos - centroid
    norm_vec = np.linalg.norm(vec)
    if norm_vec < 1e-6:
        return 90.0
    cos_angle = abs(np.dot(vec, ring_normal)) / norm_vec
    return np.degrees(np.arccos(np.clip(cos_angle, 0, 1)))


def analyze_pdb(pdb_path: str, pdb_id: str) -> List[Dict]:
    """Analyze a PDB file for sugar-protein contacts.

    Returns one dict per binding site (ligand residue).
    """
    if not HAS_BIOPYTHON:
        return []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]

    ligand_ids = LIGAND_COMP_IDS.get(pdb_id.lower(), [])
    if not ligand_ids:
        print(f"  No known ligand IDs for {pdb_id}")
        return []

    # Collect all atoms for neighbor search
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)

    # Find ligand residues
    ligand_residues = []
    protein_residues = []
    water_residues = []
    metal_atoms = []

    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            hetflag = residue.get_id()[0]
            if resname in ligand_ids:
                ligand_residues.append(residue)
            elif hetflag == " " or hetflag == "A":  # standard residues
                protein_residues.append(residue)
            elif resname == "HOH":
                water_residues.append(residue)
            elif any(a.element.strip().upper() in METAL_ELEMENTS for a in residue.get_atoms()):
                for a in residue.get_atoms():
                    if a.element.strip().upper() in METAL_ELEMENTS:
                        metal_atoms.append(a)

    results = []
    seen_sites = set()

    for lig_res in ligand_residues:
        site_key = (lig_res.get_parent().id, lig_res.get_id())
        if site_key in seen_sites:
            continue
        seen_sites.add(site_key)

        lig_atoms = list(lig_res.get_atoms())
        lig_polar = [a for a in lig_atoms if _is_polar(a)]
        lig_carbon = [a for a in lig_atoms if a.element.strip().upper() == "C"]

        # H-bonds: sugar polar ↔ protein polar within cutoff
        hbonds = []
        for la in lig_polar:
            neighbors = ns.search(la.get_vector().get_array(), HB_CUTOFF)
            for nb in neighbors:
                parent = nb.get_parent()
                if parent == lig_res:
                    continue
                if parent.get_id()[0] == " " and _is_polar(nb):
                    d = _distance(la, nb)
                    if d < HB_CUTOFF:
                        hbonds.append({
                            "sugar_atom": la.get_name(),
                            "protein_residue": f"{parent.get_resname()}{parent.get_id()[1]}",
                            "protein_atom": nb.get_name(),
                            "distance": round(d, 2),
                        })

        # Aromatic/CH-π contacts: sugar C within distance range of aromatic ring centroid
        # AND within angular cone of ring normal (geometry check)
        aromatic_contacts = []
        ch_pi_contacts = []  # geometry-verified subset
        for chain in model:
            for res in chain:
                resname = res.get_resname().strip()
                if resname not in AROMATIC_RESIDUES:
                    continue
                ring_names = AROMATIC_ATOMS.get(resname, [])
                ring6_names = AROMATIC_6RING.get(resname, ring_names)
                centroid = _ring_centroid(res, ring_names)
                normal = _ring_normal(res, ring6_names)
                if centroid is None:
                    continue
                for lc in lig_carbon:
                    lc_pos = lc.get_vector().get_array()
                    d = np.linalg.norm(lc_pos - centroid)
                    if d < CH_PI_DIST_MAX:
                        res_label = f"{resname}{res.get_id()[1]}"
                        contact = {
                            "sugar_atom": lc.get_name(),
                            "aromatic_residue": res_label,
                            "distance_to_centroid": round(d, 2),
                        }
                        aromatic_contacts.append(contact)
                        # Geometry check: is the C above/below the ring?
                        if normal is not None and d >= CH_PI_DIST_MIN:
                            angle = _ch_pi_angle(lc_pos, centroid, normal)
                            contact["angle_from_normal"] = round(angle, 1)
                            if angle < CH_PI_ANGLE_MAX:
                                ch_pi_contacts.append(contact)

        # Water bridges: water within cutoff of BOTH sugar and protein
        water_bridges = []
        for wr in water_residues:
            w_atoms = list(wr.get_atoms())
            if not w_atoms:
                continue
            wa = w_atoms[0]  # water O
            near_sugar = any(_distance(wa, la) < WATER_CUTOFF for la in lig_polar)
            if near_sugar:
                near_protein = False
                for nb in ns.search(wa.get_vector().get_array(), WATER_CUTOFF):
                    parent = nb.get_parent()
                    if parent.get_id()[0] == " " and _is_polar(nb):
                        near_protein = True
                        break
                if near_protein:
                    water_bridges.append({
                        "water_id": str(wr.get_id()[1]),
                    })

        # Metal coordination: metal within cutoff of sugar O
        metal_contacts = []
        for ma in metal_atoms:
            for la in lig_polar:
                d = _distance(ma, la)
                if d < METAL_CUTOFF:
                    metal_contacts.append({
                        "metal": ma.element.strip(),
                        "sugar_atom": la.get_name(),
                        "distance": round(d, 2),
                    })

        # Deduplicate aromatic contacts by residue
        unique_aromatic_res = set()
        for ac in aromatic_contacts:
            unique_aromatic_res.add(ac["aromatic_residue"])
        unique_chpi_res = set()
        for cc in ch_pi_contacts:
            unique_chpi_res.add(cc["aromatic_residue"])

        results.append({
            "pdb_id": pdb_id.upper(),
            "chain": lig_res.get_parent().id,
            "ligand": f"{lig_res.get_resname()}{lig_res.get_id()[1]}",
            "n_hbonds": len(hbonds),
            "n_aromatic_proximity": len(unique_aromatic_res),
            "n_ch_pi_geometry": len(unique_chpi_res),
            "n_water_bridges": len(water_bridges),
            "n_metal_contacts": len(metal_contacts),
            "aromatic_residues": sorted(unique_aromatic_res),
            "ch_pi_residues": sorted(unique_chpi_res),
            "hbond_details": hbonds[:10],
            "aromatic_details": aromatic_contacts,
            "ch_pi_details": ch_pi_contacts,
            "metal_details": metal_contacts,
        })

    return results


def analyze_all_pdbs(pdb_dir: str = "data/pdb_files") -> Dict:
    """Analyze all target PDB files."""
    pdb_path = Path(pdb_dir)
    all_results = {}

    for pdb_id in LIGAND_COMP_IDS:
        fpath = pdb_path / f"{pdb_id}.pdb"
        if not fpath.exists():
            print(f"  {pdb_id}.pdb not found — skipping")
            continue
        print(f"  Analyzing {pdb_id.upper()}...")
        results = analyze_pdb(str(fpath), pdb_id)
        all_results[pdb_id] = results

        # Print summary for first binding site
        for r in results:
            print(f"    Site {r['ligand']}: {r['n_hbonds']} HB, "
                  f"{r['n_ch_pi_geometry']}/{r['n_aromatic_proximity']} CH-π (geom/prox), "
                  f"{r['n_water_bridges']} water, "
                  f"{r['n_metal_contacts']} metal")
            if r['ch_pi_residues']:
                print(f"      CH-π verified: {', '.join(r['ch_pi_residues'])}")
            elif r['aromatic_residues']:
                print(f"      Aromatic proximity only: {', '.join(r['aromatic_residues'])}")

    return all_results


def save_contact_maps(all_results: Dict, outdir: str = "data/contact_maps"):
    """Save contact analysis to CSV and JSON."""
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # JSON with full details
    with open(f"{outdir}/plip_contacts.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"  Saved: {outdir}/plip_contacts.json")

    # CSV summary (one row per binding site)
    with open(f"{outdir}/contact_maps_auto.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "chain", "ligand", "lectin", "sugar",
                     "n_HB", "n_aromatic_proximity", "n_ch_pi_geometry",
                     "n_water_bridges", "n_metal",
                     "aromatic_residues", "ch_pi_residues"])
        for pdb_id, sites in sorted(all_results.items()):
            info = {
                "5cna": ("ConA", "MeαMan"), "1cvn": ("ConA", "trimannose"),
                "2uvo": ("WGA", "GlcNAc"), "2pel": ("PNA", "galactose"),
                "1a3k": ("Galectin-3", "LacNAc"), "1lzb": ("Lysozyme", "chitobiose"),
            }
            lectin, sugar = info.get(pdb_id, ("?", "?"))
            for site in sites:
                w.writerow([
                    site["pdb_id"], site["chain"], site["ligand"],
                    lectin, sugar,
                    site["n_hbonds"], site["n_aromatic_proximity"],
                    site["n_ch_pi_geometry"],
                    site["n_water_bridges"], site["n_metal_contacts"],
                    "; ".join(site["aromatic_residues"]),
                    "; ".join(site["ch_pi_residues"]),
                ])
    print(f"  Saved: {outdir}/contact_maps_auto.csv")


if __name__ == "__main__":
    print("Contact Analysis — BioPython-based PDB analyzer")
    print("=" * 60)
    results = analyze_all_pdbs()
    if results:
        save_contact_maps(results)
    else:
        print("No results. Ensure PDB files are in data/pdb_files/")