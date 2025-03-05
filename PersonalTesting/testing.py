from molecule import Molecule
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw

testAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
testBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]
mol1 = Molecule()
mol1.fromAtomsAndBonds(testAtoms, testBonds)

print("Graph Nodes (Atoms):")
print(mol1.graph.nodes(data=True)) 
print("\nGraph Edges (Bonds):")
print(mol1.graph.edges(data=True))  

sdfFile = "/home/codingnewt/miniconda3/chem274a/final-assignment-ProgrammingNewt/sdf_files/sdf/AAKJLRGGTJKAMG-UHFFFAOYSA-N.sdf"
mol2 = Molecule()
try:
    mol2.fromSdf(sdfFile, includeHydrogen=False)
    print("\nSDF Graph Nodes (Atoms):")
    print(mol2.graph.nodes(data=True)) 
    print("\nSDF Graph Edges (Bonds):")
    print(mol2.graph.edges(data=True))
except FileNotFoundError:
    print(f"File {sdfFile} not found.")
except ValueError as e:
    print(f"Error parsing SDF file, {e}")
try:
    fingerprint1 = mol1.generateFingerprint(pathLength=7, fingerprintSize=2048, bitsPerHash=2)
    print("\nGenerated Fingerprint for mol1 (First 50 bits):")
    print(fingerprint1[:50])

    fingerprint2 = mol2.generateFingerprint(pathLength=7, fingerprintSize=2048, bitsPerHash=2)
    print("\nGenerated Fingerprint for mol2 (First 50 bits):")
    print(fingerprint2[:50])
except Exception as e:
    print(f"Error generating fingerprint: {e}")

try:
    aromaticRings1 = mol1.detectAromaticity()
    print("\nDetected Aromatic Rings in mol1:")
    for ring in aromaticRings1:
        print(f"Aromatic Ring: {ring}")

    aromaticRings2 = mol2.detectAromaticity()
    print("\nDetected Aromatic Rings in mol2:")
    for ring in aromaticRings2:
        print(f"Aromatic Ring: {ring}")
except Exception as e:
    print(f"Error detecting aromaticity: {e}")
try:
    mol3 = Molecule()
    mol3.fromAtomsAndBonds(testAtoms, testBonds)
    print("\nComparing mol1 and mol3 (should be equal):", mol1 == mol3)
except Exception as e:
    print(f"Error testing equality operator: {e}")
try:
    print("\nComparing mol1 and mol2 (should not be equal):", mol1 == mol2)
except Exception as e:
    print(f"Error testing equality operator: {e}")
rdkitMol = Chem.SDMolSupplier(sdfFile, removeHs=False)[0]
if rdkitMol is None:
    print(f"Failed to load molecule from {sdfFile}")
else:
    aromaticAtoms = []
    for atom in rdkitMol.GetAtoms():
        if atom.GetIsAromatic():
            aromaticAtoms.append(atom.GetIdx())

    img = Draw.MolToImage(rdkitMol, highlightAtoms=aromaticAtoms)
    img.save("molecule_visualization.png")
    print("Molecule visualization saved as 'molecule_visualization.png'")
try:
    subAtoms = {'C1': 'C', 'C2': 'C'}
    subBonds = [('C1', 'C2', 1)]
    substructure = Molecule()
    substructure.fromAtomsAndBonds(subAtoms, subBonds)

    print("\nDoes mol2 contain the substructure? (should be True):", mol2.containsSubstructure(substructure))
    print("\nDoes mol1 contain the substructure? (should be True):", mol1.containsSubstructure(substructure))
except Exception as e:
    print(f"Error testing substructure screening: {e}")

