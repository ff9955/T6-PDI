import numpy as np

# tolerance for floating-point comparison
TOL = 1e-6


def read_xyz(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    # Keep header lines
    header = lines[:2]
    atoms = []
    for i, line in enumerate(lines[2:]):
        parts = line.split()
        if len(parts) < 4:
            continue
        atom = parts[0]
        coords = list(map(float, parts[1:4]))
        atoms.append((atom, np.array(coords), i + 1, line))  # store 1-based index + raw line
    return header, atoms


def match_atoms(file1, file2, output_file):
    header1, atoms1 = read_xyz(file1)
    header2, atoms2 = read_xyz(file2)

    matches = []
    used2 = set()

    for atom1, coords1, idx1, _ in atoms1:
        found_match = False
        for j, (atom2, coords2, idx2, _) in enumerate(atoms2):
            if j in used2:
                continue
            if atom1 == atom2 and np.linalg.norm(coords1 - coords2) < TOL:
                matches.append((idx1, idx2))
                used2.add(j)
                found_match = True
                break
        if not found_match:
            print(f"Warning: no match found for atom at line {idx1} in {file1}")

    with open(output_file, "w") as out:
        for i1, i2 in matches:
            out.write(f"{i1} {i2}\n")

    return matches, header1, atoms1, atoms2


def reorder_xyz_and_check(file1, file2, mapping_file, output_file):
    matches, header1, atoms1, atoms2 = match_atoms(file1, file2, mapping_file)

    # Build lookup for file2 atoms by index
    atom2_dict = {idx2: line for _, _, idx2, line in atoms2}

    # Reorder according to file1 order
    reordered_lines = []
    for idx1, idx2 in matches:
        reordered_lines.append(atom2_dict[idx2])

    with open(output_file, "w") as out:
        out.writelines(header1)
        out.writelines(reordered_lines)

def reorder_xyz_new(mol_file, mapping_file):

    new_coords = []

    full_xyz = open(mol_file).readlines()
    header = full_xyz[:2]
    mol_coords = full_xyz[2:]

    map_elements = np.loadtxt(mapping_file)

    for element in map_elements:
        new_coords.append(mol_coords[int(element[1])-1])

    with open(f'{mol_file[:-4]}-mapped.xyz', 'w') as out:
        out.writelines(header)
        out.writelines(new_coords)
    

if __name__ == "__main__":
    file1 = "mol4-mapped.xyz"
    file2 = "mol4.xyz"
    mapping_file = "mols-14-mapping.txt"
    reordered_file = "mol4-map-test.xyz"

    reorder_xyz_and_check(file1, file2, mapping_file, reordered_file)
    print(f"Mapping written to {mapping_file}")
    print(f"Reordered XYZ written to {reordered_file}")

#    reorder_xyz_new('mol4.xyz', 'mols-13-mapping.txt')


