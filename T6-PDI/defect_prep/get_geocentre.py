import numpy as np

# Define atomic masses for common elements (in atomic mass units)
atomic_masses = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    # Add more as needed
}

def read_xyz(filename):
    """Reads an XYZ file and returns a list of tuples: (element, x, y, z)"""
    atoms = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[0].strip())
        for line in lines[2:2+num_atoms]:
            parts = line.split()
            element = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((element, x, y, z))
    return atoms

def center_of_mass(atoms, selection=None):
    """
    Computes the center of mass for a list of atoms.
    atoms: list of tuples (element, x, y, z)
    selection: list of indices (0-based) to include; if None, use all atoms
    """
    if selection is None:
        selection = range(len(atoms))

    total_mass = 0.0
    com = np.array([0.0, 0.0, 0.0])

    for i in selection:
        element, x, y, z = atoms[i]
        mass = atomic_masses.get(element)
        if mass is None:
            raise ValueError(f"Mass for element '{element}' not defined.")
        total_mass += mass
        com += mass * np.array([x, y, z])

    com /= total_mass
    #convert to Bohr radii
    com = com*1.8897259886
    return com

# Example usage
if __name__ == "__main__":
    xyz_file = "final_nvt_geom_defect2.xyz"  #Replace with your file
    atoms = read_xyz(xyz_file)

    # Optional: select only certain atoms, e.g., [0, 1, 2]
    selection = [13402, 13403, 13404, 13414, 13412, 13413]

    com = center_of_mass(atoms, selection)
    print(f"Center of mass (x, y, z): {com[0]:.6f}, {com[1]:.6f}, {com[2]:.6f}")

