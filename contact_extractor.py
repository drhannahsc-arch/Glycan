"""
Contact extractor for lectin-glycan PDB structures.

Takes a PDB file and identifies all protein-sugar contacts:
  - H-bonds (O/N...O/N < 3.5 Å)
  - CH-π contacts (sugar C near aromatic ring atoms < 4.5 Å)
  - Metal coordination (sugar O near Ca²⁺/Mn²⁺ < 3.0 Å)
  - Bridging waters (water O within 3.5 Å of BOTH sugar and protein polar atoms)
  - Buried hydroxyl classification (equatorial vs axial)

Requires: biopython (pip install biopython)

Data source: PDB files from EBI (www.ebi.ac.uk/pdbe/entry-files/)
"""

import os
import urllib.request
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple
from pathlib import Path

from Bio.PDB import PDBParser

# =====================================================================
# ContactMap dataclass — the interface between extraction and scoring
# =====================================================================

@dataclass
class ContactMap:
    """All contacts between a sugar ligand and its binding site.

    This is THE interface between the contact extractor and the scorer.
    Students fill this from PDB parsing; the scorer consumes it.
    """
    pdb_id: str
    chain: str
    sugar_resname: str
    sugar_resid: int

    # Contact counts
    n_hbonds: int = 0              # Sugar O/N ... Protein O/N < 3.5 Å
    n_ch_pi: int = 0              # Sugar C ... Aromatic ring atoms < 4.5 Å
    n_buried_eq_OH: int = 0       # Equatorial OHs in contact with protein
    n_buried_ax_OH: int = 0       # Axial OHs in contact with protein
    n_buried_NAc: int = 0         # NAc groups in contact with protein
    n_frozen_torsions: int = 0    # Glycosidic linkages frozen (0 for mono)
    n_water_bridges: int = 0      # Waters bridging sugar and protein

    # Residue-level detail
    aromatic_residues: List[str] = field(default_factory=list)
    metal_ions: List[str] = field(default_factory=list)
    hbond_partners: List[str] = field(default_factory=list)

    # Raw contact list for debugging
    contacts: List[dict] = field(default_factory=list)

    def summary(self) -> str:
        """One-line summary."""
        return (
            f"{self.pdb_id} chain {self.chain} {self.sugar_resname}#{self.sugar_resid}: "
            f"HB={self.n_hbonds} CH-π={self.n_ch_pi} "
            f"eq_OH={self.n_buried_eq_OH} ax_OH={self.n_buried_ax_OH} "
            f"NAc={self.n_buried_NAc} waters={self.n_water_bridges} "
            f"metals={len(self.metal_ions)} "
            f"aromatics={self.aromatic_residues}"
        )


# =====================================================================
# Sugar residue classification
# =====================================================================

# PDB 3-letter codes for common sugars
SUGAR_RESNAMES = {
    # Mannose
    "MAN", "BMA", "MMA",   # alpha-Man, beta-Man, methyl-alpha-Man
    # Glucose
    "GLC", "BGC", "AGC",   # alpha-Glc, beta-Glc, ...
    # Galactose
    "GAL", "GLA",           # beta-Gal, alpha-Gal
    # GlcNAc
    "NAG", "NDG",           # N-acetyl-D-glucosamine variants
    # GalNAc
    "NGA", "A2G",           # N-acetyl-D-galactosamine
    # Fucose
    "FUC", "FCA",
    # Sialic acid
    "SIA", "SLB",
    # Xylose
    "XYS", "XYP",
    # Generic
    "TRE", "LAT", "SUC", "MAL",  # disaccharides
}

# Sugar type classification for scoring
SUGAR_TYPE = {
    "MAN": "Man", "BMA": "Man", "MMA": "Man",
    "GLC": "Glc", "BGC": "Glc", "AGC": "Glc",
    "GAL": "Gal", "GLA": "Gal",
    "NAG": "GlcNAc", "NDG": "GlcNAc",
    "NGA": "GalNAc", "A2G": "GalNAc",
    "FUC": "Fuc", "FCA": "Fuc",
    "SIA": "NeuAc", "SLB": "NeuAc",
}

# Aromatic residues (for CH-π detection)
AROMATIC_RESIDUES = {"TRP", "TYR", "PHE", "HIS"}

# Aromatic ring atom names per residue
AROMATIC_RING_ATOMS = {
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}

# Metal ions relevant to lectin binding
METAL_RESNAMES = {"CA", "MN", "ZN", "MG", "FE"}

# Hydroxyl classification by sugar type and atom name
# These map sugar OH oxygen atoms to equatorial/axial based on standard pyranose geometry
# Convention: for 4C1 chair (the dominant conformation)
HYDROXYL_CLASS = {
    # Mannose: C2-OH is axial, rest equatorial
    "Man": {"O2": "ax", "O3": "eq", "O4": "eq", "O6": "eq"},
    # Glucose: all equatorial
    "Glc": {"O2": "eq", "O3": "eq", "O4": "eq", "O6": "eq"},
    # Galactose: C4-OH is axial, rest equatorial
    "Gal": {"O2": "eq", "O3": "eq", "O4": "ax", "O6": "eq"},
    # GlcNAc: like Glc but C2 has NAc instead of OH
    "GlcNAc": {"O3": "eq", "O4": "eq", "O6": "eq", "N2": "NAc"},
    # GalNAc: like Gal but C2 has NAc
    "GalNAc": {"O3": "eq", "O4": "ax", "O6": "eq", "N2": "NAc"},
    # Fucose: C4-OH axial (6-deoxy galactose configuration)
    "Fuc": {"O2": "eq", "O3": "eq", "O4": "ax"},
}


# =====================================================================
# PDB fetching
# =====================================================================

PDB_CACHE_DIR = Path(__file__).parent / "pdb_cache"

def fetch_pdb(pdb_id: str, cache_dir: Optional[str] = None) -> str:
    """Download PDB file from EBI if not cached. Returns file path."""
    cache = Path(cache_dir) if cache_dir else PDB_CACHE_DIR
    cache.mkdir(exist_ok=True)

    pdb_lower = pdb_id.lower()
    fn = cache / f"pdb{pdb_lower}.ent"

    if fn.exists():
        return str(fn)

    url = f"https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_lower}.ent"
    try:
        data = urllib.request.urlopen(url, timeout=30).read()
        with open(fn, "wb") as f:
            f.write(data)
        return str(fn)
    except Exception as e:
        raise RuntimeError(f"Cannot fetch PDB {pdb_id} from EBI: {e}")


# =====================================================================
# Contact extraction
# =====================================================================

# Distance cutoffs (Angstroms)
HB_CUTOFF = 3.2       # Heavy-atom donor...acceptor (strong H-bonds only)
CHPI_CUTOFF = 4.5     # Sugar C to aromatic ring atom
METAL_CUTOFF = 3.0    # Sugar O to metal ion
WATER_BRIDGE_CUTOFF = 3.5  # Water O to polar atom (both sugar and protein)
CONTACT_CUTOFF = 4.0  # General contact cutoff


def extract_contacts(pdb_id: str, sugar_resname: Optional[str] = None,
                     chain: Optional[str] = None,
                     pdb_file: Optional[str] = None,
                     cache_dir: Optional[str] = None) -> List[ContactMap]:
    """Extract all protein-sugar contacts from a PDB structure.

    Args:
        pdb_id: 4-character PDB code (e.g., "3ZSJ")
        sugar_resname: Filter to specific sugar (e.g., "GAL"). None = all sugars.
        chain: Filter to specific chain. None = all chains.
        pdb_file: Path to local PDB file. If None, fetches from EBI.
        cache_dir: Cache directory for downloaded PDB files.

    Returns:
        List of ContactMap objects, one per sugar residue found.
    """
    # Load structure
    if pdb_file is None:
        pdb_file = fetch_pdb(pdb_id, cache_dir)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]

    # Collect atoms by category
    protein_atoms = []
    sugar_residues = []
    water_atoms = []
    metal_atoms = []

    for ch in model:
        for res in ch.get_residues():
            resname = res.resname.strip()
            het_flag = res.id[0]

            if het_flag == " ":
                # Standard protein residue
                for a in res.get_atoms():
                    protein_atoms.append(a)
            elif resname == "HOH":
                for a in res.get_atoms():
                    water_atoms.append(a)
            elif resname in METAL_RESNAMES:
                for a in res.get_atoms():
                    metal_atoms.append(a)
            elif resname in SUGAR_RESNAMES:
                if sugar_resname and resname != sugar_resname:
                    continue
                if chain and ch.id != chain:
                    continue
                sugar_residues.append(res)

    if not sugar_residues:
        # Try broader matching — some PDB files use non-standard codes
        for ch in model:
            for res in ch.get_residues():
                if res.id[0] not in (" ", "W") and res.resname.strip() not in METAL_RESNAMES:
                    if res.resname.strip() != "HOH":
                        sugar_residues.append(res)

    results = []

    for sugar_res in sugar_residues:
        resname = sugar_res.resname.strip()
        resid = sugar_res.id[1]
        sugar_chain = sugar_res.get_parent().id
        sugar_type = SUGAR_TYPE.get(resname, "unknown")

        cm = ContactMap(
            pdb_id=pdb_id,
            chain=sugar_chain,
            sugar_resname=resname,
            sugar_resid=resid,
        )

        sugar_atoms = list(sugar_res.get_atoms())
        sugar_polar = [a for a in sugar_atoms if a.element in ("O", "N")]
        sugar_carbon = [a for a in sugar_atoms if a.element == "C"]

        # ── H-bonds: sugar O/N ... protein O/N < 3.5 Å ──────────
        hbond_pairs = set()
        for sa in sugar_polar:
            for pa in protein_atoms:
                if pa.element not in ("O", "N"):
                    continue
                d = sa - pa
                if d < HB_CUTOFF:
                    pres = pa.get_parent()
                    pair_key = (sa.name, pres.resname, pres.id[1], pa.name)
                    if pair_key not in hbond_pairs:
                        hbond_pairs.add(pair_key)
                        cm.hbond_partners.append(
                            f"{sa.name}...{pres.resname}{pres.id[1]}.{pa.name} ({d:.1f}A)"
                        )
                        cm.contacts.append({
                            "type": "hbond",
                            "sugar_atom": sa.name,
                            "protein_res": f"{pres.resname}{pres.id[1]}",
                            "protein_atom": pa.name,
                            "distance_A": round(d, 2),
                        })
        cm.n_hbonds = len(hbond_pairs)

        # ── CH-π: sugar C ... aromatic ring atoms < 4.5 Å ───────
        ch_pi_residues = set()
        ch_pi_count = 0
        for sc in sugar_carbon:
            for pa in protein_atoms:
                pres = pa.get_parent()
                if pres.resname not in AROMATIC_RESIDUES:
                    continue
                ring_atoms = AROMATIC_RING_ATOMS.get(pres.resname, set())
                if pa.name not in ring_atoms:
                    continue
                d = sc - pa
                if d < CHPI_CUTOFF:
                    res_key = f"{pres.resname}{pres.id[1]}"
                    if res_key not in ch_pi_residues:
                        ch_pi_residues.add(res_key)
                        cm.aromatic_residues.append(res_key)
                    ch_pi_count += 1

        # Count effective CH-π contacts (unique aromatic residues, not atom pairs)
        cm.n_ch_pi = len(ch_pi_residues)

        # ── Metal coordination: sugar O ... metal < 3.0 Å ───────
        for sa in sugar_polar:
            for ma in metal_atoms:
                d = sa - ma
                if d < METAL_CUTOFF:
                    metal_res = ma.get_parent()
                    metal_key = f"{metal_res.resname}{metal_res.id[1]}"
                    if metal_key not in cm.metal_ions:
                        cm.metal_ions.append(metal_key)
                    cm.contacts.append({
                        "type": "metal",
                        "sugar_atom": sa.name,
                        "metal": metal_key,
                        "distance_A": round(d, 2),
                    })

        # ── Bridging waters ──────────────────────────────────────
        for wa in water_atoms:
            near_sugar = False
            near_protein = False
            for sa in sugar_polar:
                if sa - wa < WATER_BRIDGE_CUTOFF:
                    near_sugar = True
                    break
            if not near_sugar:
                continue
            for pa in protein_atoms:
                if pa.element not in ("O", "N"):
                    continue
                if pa - wa < WATER_BRIDGE_CUTOFF:
                    near_protein = True
                    break
            if near_sugar and near_protein:
                cm.n_water_bridges += 1

        # ── Classify buried hydroxyls ────────────────────────────
        oh_class = HYDROXYL_CLASS.get(sugar_type, {})
        for sa in sugar_polar:
            aname = sa.name.strip()
            if aname not in oh_class:
                continue
            # Check if this OH is in contact with protein
            buried = False
            for pa in protein_atoms:
                if sa - pa < CONTACT_CUTOFF:
                    buried = True
                    break
            if buried:
                cls = oh_class[aname]
                if cls == "eq":
                    cm.n_buried_eq_OH += 1
                elif cls == "ax":
                    cm.n_buried_ax_OH += 1
                elif cls == "NAc":
                    cm.n_buried_NAc += 1

        results.append(cm)

    return results


def extract_contacts_simple(pdb_id: str, sugar_resname: str,
                            chain: Optional[str] = None) -> ContactMap:
    """Convenience wrapper: extract contacts for a single sugar.

    Returns the first matching ContactMap, or raises ValueError.
    """
    results = extract_contacts(pdb_id, sugar_resname=sugar_resname, chain=chain)
    if not results:
        raise ValueError(f"No {sugar_resname} found in {pdb_id}" +
                         (f" chain {chain}" if chain else ""))
    return results[0]


# =====================================================================
# Batch extraction for reference structures
# =====================================================================

REFERENCE_EXTRACTIONS = [
    # (pdb_id, sugar_resname, chain, description)
    ("3ZSJ", "GAL", None, "Galectin-3 + Gal (from LacNAc)"),
    ("5CNA", "MMA", None, "ConA + methyl-alpha-mannopyranoside"),
    ("1SL5", "GAL", None, "DC-SIGN + sugar"),
    ("2UVO", "NAG", None, "WGA + GlcNAc"),
]


def extract_reference_set(cache_dir: Optional[str] = None) -> Dict[str, ContactMap]:
    """Extract contacts for all reference lectin structures.

    Returns dict keyed by "PDBID_RESNAME".
    """
    results = {}
    for pdb_id, resname, chain, desc in REFERENCE_EXTRACTIONS:
        try:
            cms = extract_contacts(pdb_id, sugar_resname=resname,
                                   chain=chain, cache_dir=cache_dir)
            if cms:
                key = f"{pdb_id}_{resname}"
                results[key] = cms[0]
        except Exception as e:
            print(f"WARNING: {pdb_id} {resname}: {e}")
    return results


if __name__ == "__main__":
    import sys

    cache = "/home/claude/pdb_cache"

    print("Contact Extractor — Reference Structure Analysis")
    print("=" * 65)

    test_cases = [
        ("3ZSJ", "GAL", None),
        ("5CNA", "MMA", None),
        ("1SL5", None, None),
        ("2UVO", "NAG", None),
    ]

    for pdb_id, resname, chain in test_cases:
        print(f"\n--- {pdb_id} {resname or 'ALL'} ---")
        try:
            cms = extract_contacts(pdb_id, sugar_resname=resname,
                                   chain=chain, cache_dir=cache)
            for cm in cms[:3]:  # limit output
                print(f"  {cm.summary()}")
                if cm.hbond_partners:
                    print(f"    H-bonds: {', '.join(cm.hbond_partners[:5])}")
        except Exception as e:
            print(f"  ERROR: {e}")