import os
import urllib.request
from Bio.PDB import PDBParser

instance_id = "1GPV_A"

pdb_id = instance_id.split("_")[0]
chain_id = instance_id.split("_")[1]

# Step 1: Download the 1GPV PDB file
filename = os.path.join("instances", f"{pdb_id}.pdb")
if not os.path.exists(filename):
    urllib.request.urlretrieve(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", filename
    )


# Step 2: Parse the PDB file and extract chain A atoms
parser = PDBParser()
structure = parser.get_structure(pdb_id, filename)
model = structure[0]
chain = model[chain_id]

atoms = []
for residue in chain:
    res_seq = residue.id[1]  # Get the residue sequence identifier
    res_name = residue.resname
    for atom in residue:
        atoms.append(
            {
                "serial": atom.serial_number,
                "atom_name": atom.name.strip(),
                "residue_name": res_name,
                "res_seq": res_seq,
                "x": atom.coord[0],
                "y": atom.coord[1],
                "z": atom.coord[2],
            }
        )

# Step 3: Sort the atoms by residue sequence identifier and atom serial number
sorted_atoms = sorted(atoms, key=lambda a: (a["res_seq"], a["serial"]))

# Step 4: Save the sorted atoms to 1GPV_A.txt
filename = os.path.join("instances", f"{instance_id}.txt")
with open(filename, "w") as f:
    f.write(
        "ATOM_SERIAL_NUMBER,ATOM_NAME,RESIDUE_NAME,RESIDUE_SEQUENCE_IDENTIFIER,CHAIN_ID,X,Y,Z\n"
    )
    for atom in sorted_atoms:
        line = f"{atom['serial']},{atom['atom_name']},{atom['residue_name']},{atom['res_seq']},A,{atom['x']:.8f},{atom['y']:.8f},{atom['z']:.8f}\n"
        f.write(line)

# Steps 5-7: Generate instance files for each N and epsilon
Ns = [100, 200]
epsilons = [0.04, 0.08, 0.12, 0.16]

for N in Ns:
    V = sorted_atoms[:N]
    E = []
    # Generate instance solution
    filename = os.path.join("instances", f"{instance_id}_N_{N}.sol")
    with open(filename, "w") as f:
        f.write(f"ATOM_SERIAL_NUMBER,X,Y,Z\n")
        for i in range(len(V)):
            atom_i = V[i]
            f.write(
                f'{atom_i["serial"]} {atom_i["x"]:.8f} {atom_i["y"]:.8f} {atom_i["z"]:.8f}\n'
            )

    # Generate pairs (i, j) where i and j are in same or consecutive residues
    for i in range(len(V)):
        atom_i = V[i]
        res_i = atom_i["res_seq"]
        for j in range(i + 1, len(V)):
            atom_j = V[j]
            res_j = atom_j["res_seq"]
            if res_j == res_i or res_j == res_i + 1:
                # Calculate Euclidean distance
                dx = atom_i["x"] - atom_j["x"]
                dy = atom_i["y"] - atom_j["y"]
                dz = atom_i["z"] - atom_j["z"]
                distance = (dx**2 + dy**2 + dz**2) ** 0.5
                E.append((atom_i["serial"], atom_j["serial"], distance))
            elif res_j > res_i + 1:
                break  # No need to check further

    # Create instance files for each epsilon
    for epsilon in epsilons:
        filename = os.path.join(
            "instances", f"{instance_id}_N_{N}_eps_{epsilon:.2f}.txt"
        )
        with open(filename, "w") as f:
            f.write(
                "ATOM_SERIAL_NUMBER_I,ATOM_SERIAL_NUMBER_J,DISTANCE_IJ,LOWER_BOUND_IJ,UPPER_BOUND_IJ\n"
            )
            for i_serial, j_serial, d in E:
                lower = (1 - epsilon) * d
                upper = (1 + epsilon) * d
                f.write(f"{i_serial},{j_serial},{d:.8f},{lower:.8f},{upper:.8f}\n")
