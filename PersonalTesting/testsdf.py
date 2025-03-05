from provided import parse_sdf

sdfFile = "/home/codingnewt/miniconda3/chem274a/final-assignment-ProgrammingNewt/sdf_files/sdf/AAKJLRGGTJKAMG-UHFFFAOYSA-N.sdf"
try:
    atoms, bonds = parse_sdf(sdfFile, include_hydrogen=False)
    print("Atoms:")
    print(atoms)
    print("\nBonds:")
    print(bonds) 
    
except FileNotFoundError:
    print(f"File {sdfFile} not found.")
except ValueError as e:
    print(f"Error parsing {sdfFile}: {e}")
