import pytest
from molecule import Molecule
import os

def testFromAtomsAndBonds():
    testAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
    testBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]
    mol = Molecule()
    mol.fromAtomsAndBonds(testAtoms, testBonds)

    for atom in ['C1', 'C2', 'N1']:
        assert mol.graph.nodes[atom]['element'] == testAtoms[atom]

    for bond in testBonds:
        atom1, atom2, bondOrder = bond
        assert (atom1, atom2) in mol.graph.edges
        assert mol.graph[atom1][atom2]['bondOrder'] == bondOrder

def testFromSdf():
    sdfFile = "sdf_files/sdf/AAKJLRGGTJKAMG-UHFFFAOYSA-N.sdf"
    mol = Molecule()
    mol.fromSdf(sdfFile)

    assert 'C1' in mol.graph.nodes
    assert 'N1' in mol.graph.nodes
    assert ('C1', 'C2') in mol.graph.edges

def testGenerateFingerprint():
    testAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
    testBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]
    mol = Molecule()
    mol.fromAtomsAndBonds(testAtoms, testBonds)

    fingerprint = mol.generateFingerprint(pathLength=7, fingerprintSize=2048, bitsPerHash=2)
    assert len(fingerprint) == 2048

    bitCount = sum(fingerprint)
    assert bitCount > 0  

def testDetectAromaticity():
    sdfFile = "sdf_files/sdf/AAKJLRGGTJKAMG-UHFFFAOYSA-N.sdf"
    mol = Molecule()
    mol.fromSdf(sdfFile)

    aromaticRings = mol.detectAromaticity()
    for ring in aromaticRings:
        assert len(ring) >= 6  
    assert len(aromaticRings) > 0

def testEquality():
    testAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
    testBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]

    mol1 = Molecule()
    mol1.fromAtomsAndBonds(testAtoms, testBonds)

    mol2 = Molecule()
    mol2.fromAtomsAndBonds(testAtoms, testBonds)

    assert mol1 == mol2

def testSubstructureScreening():
    mainAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
    mainBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]
    subAtoms = {'C1': 'C', 'C2': 'C'}
    subBonds = [('C1', 'C2', 1)]

    mainMolecule = Molecule()
    mainMolecule.fromAtomsAndBonds(mainAtoms, mainBonds)

    substructure = Molecule()
    substructure.fromAtomsAndBonds(subAtoms, subBonds)

    result = mainMolecule.containsSubstructure(substructure)
    assert result

def testVisualization():
    
    testAtoms = {'C1': 'C', 'C2': 'C', 'N1': 'N'}
    testBonds = [('C1', 'C2', 1), ('C2', 'N1', 2)]
    mol = Molecule()
    mol.fromAtomsAndBonds(testAtoms, testBonds)

    savePath = "testMol.png"
    mol.visualize(savePath=savePath)

    assert os.path.exists(savePath)
