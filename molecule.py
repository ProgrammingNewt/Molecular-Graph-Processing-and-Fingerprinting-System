import hashlib
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from provided import parse_sdf
import matplotlib.pyplot as plt
class Molecule:
    def __init__(self):
        """
        Molecule class initalization
        """
        self.graph = nx.Graph()
        self.atoms = {}  #store atom names and elements
        self.bonds = []  # store bond information
    def fromAtomsAndBonds(self, atoms: dict, bonds: list):
            """
            Makes a "Molecule" object from the given bonds and atoms

            Parameters:
            ----------
            atoms: dict
            A dict with atom names as keys and the elements as its vals
            bonds: list
            A list of storing the bonds
            """
            self.atoms = atoms
            self.bonds = bonds
            for atomName, element in atoms.items():
                
                self.graph.add_node(atomName, element=element)

            for atom1, atom2, bondOrder in bonds:
                
                self.graph.add_edge(atom1, atom2, bondOrder=bondOrder)

    def fromSdf(self, fileName: str, includeHydrogen: bool = False):
        """
        Makes the "Molecule" object from the atom and bond data within an sdf file
        Parameters:
        ----------
        fileName: str
            path or name of sdf file
        includeHydrogen: bool
            Include hydrograns, needed for built in parse_sdf function.
        """
        atoms, bonds = parse_sdf(fileName, include_hydrogen=includeHydrogen)
        self.fromAtomsAndBonds(atoms, bonds)


    def generateFingerprint(self, pathLength=7, fingerprintSize=2048, bitsPerHash=2):
        """
        Makes a fingerprint for each molecule

        Parameters:
        ----------
        pathLength: int
            max path length
        fingerprintSize: int
            number of bits for fingerprint size.
        bitsPerHash: int
            bits for each hash

        Returns:
        -------
        list
            A fingerprint vector
        """
        fingerprint = [0] * fingerprintSize

        for startNode in self.graph.nodes:
            paths = generatePaths(self.graph, startNode, pathLength)
            for path in paths:
                pathStr = ""
                for atom in path:
                    pathStr += f"{atom}:{self.graph.nodes[atom]['element']}-"
                pathStr = pathStr.rstrip("-")

                hashValue = int(hashlib.sha256(pathStr.encode()).hexdigest(), 16)
                for i in range(bitsPerHash):
                    index = (hashValue + i) % fingerprintSize
                    fingerprint[index] = 1

        return fingerprint
    
    def detectAromaticity(self):
        """
        find aromatic rings in struct
        
        Returns:
        -------
        list
            A list of aromatic rings
        """
        aromaticRings = []
        rings = nx.cycle_basis(self.graph)  # Find all rings in the graph

        for ring in rings:
            if self.isRingAromatic(ring):
                aromaticRings.append(ring)
                for atom in ring:
                    self.graph.nodes[atom]['aromatic'] = True

        return aromaticRings
    
    def isRingAromatic(self, ring):
        if (len(ring) < 6):
            return False
        for atom in ring:
            if (self.graph.nodes[atom]['element'] not in ['C', 'N', 'O', 'S']):
                return False
        piElectrons = 0
        
        for i in range(len(ring)):
            atom = ring[i]
            nextAtom = ring[(i + 1) % len(ring)]  # Wrap to form a ring
            bondOrder = self.graph[atom][nextAtom]['bondOrder']
            if (bondOrder == 2):
                piElectrons += 2
            elif (self.graph.nodes[atom]['element'] in ['N', 'O', 'S']):
                totalBonds = sum(self.graph[atom][neighbor]['bondOrder'] for neighbor in self.graph.neighbors(atom))
                if (totalBonds <= 2): 
                    piElectrons += 2
                    
        if ((piElectrons - 2) % 4 == 0):
            return True
        return False
    
    
    def __eq__(self, mol2):
        """
        sees if two mols are the same based on their unique fingerprints

        Parameters:
        ----------
        other: Molecule
            The other molecule object

        Returns:
        -------
        bool
            True or falsed if they're the same
        """
        mol1hash = self.generateFingerprint()
        mol2hash = mol2.generateFingerprint()

        return mol1hash == mol2hash
    
    
    
    
    def containsSubstructure(self, substructure):
        """
        Sees if the substrcuture exists within larger structure

        Parameters:
        ----------
        substructure: Molecule
            The substructure we're checking for

        Returns:
        -------
        bool
            True if the substructure is there, false if not
        """
        def nMatch(a, b):
            return a['element'] == b['element']
        def eMatches(a, b):
            return a['bondOrder'] == b['bondOrder']
        
        doesMatch = GraphMatcher(self.graph, substructure.graph, node_match=nMatch, edge_match=eMatches)
        matchStruct = doesMatch.subgraph_is_isomorphic()

        return matchStruct
    
    
    def visualize(self, savePath=None):
       
        cpkColors = {
            'H': 'white',
            'C': 'green',
            'N': 'blue',
            'O': 'red',
            'S': 'yellow'
        }
        otherAtoms = 'gray'
        nodeColors = [
            cpkColors.get(self.graph.nodes[node]['element'], otherAtoms)
            for node in self.graph.nodes]

        plt.figure(figsize=(6, 6))
        pos = nx.spring_layout(self.graph)

        nx.draw(
            self.graph,pos,with_labels=True,node_color=nodeColors,edge_color='purple',node_size=600,font_size=8
        )

        if savePath:
            plt.savefig(savePath)
        else:
            plt.show()


def generatePaths(graph, startNode, maxLength):
    """
    creates all paths from one node to another for a given len

    Parameters:
    ----------
    graph: 
        The molecule graph object
    startNode: str
        The starting node
    maxLength: int
        The max path len

    Returns:
    -------
    list
        A list of paths
    """
    paths = []
    def explorePath(currentPath, remainingLength):
        #recursive DFS
        if (remainingLength > 0):
            currentNode = currentPath[-1]
            neighbors = graph.neighbors(currentNode)
            
            for neighbor in neighbors:
                
                if (neighbor not in currentPath):
                    newPath = currentPath + [neighbor]
                    paths.append(newPath)
                    explorePath(newPath, remainingLength - 1)
                    
    explorePath([startNode], maxLength)
    return paths
