import sys
from molecule import Molecule

mainMolecule = Molecule()
mainMolecule.fromSdf(sys.argv[1])
substructure = Molecule()
substructure.fromSdf(sys.argv[2])

result = mainMolecule.containsSubstructure(substructure)
if result:
    print("Substructure was found")
else:
    print("Substructure was not found :(")
