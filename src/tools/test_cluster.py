from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import torch

def load_molecules():
    routes = []
    route_file_name = '../dataset/retro190.txt'
    for line in open(route_file_name, "r"):
        routes.append(line)
    print('%d routes extracted from %s loaded' % (len(routes), route_file_name))
    return routes

def preprocess(X, fp_dim):
    # Compute fingerprint from mol to feature
    mol = Chem.MolFromSmiles(X)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=int(fp_dim),useChirality=True)
    onbits = list(fp.GetOnBits())
    arr = np.zeros(fp.GetNumBits())
    arr[onbits] = 1
    return arr

target_products = load_molecules()
target_products_fps = [preprocess(mol, 256) for mol in target_products]


import faiss                     # make faiss available
## to numpy
target_products_fps_np = np.array(target_products_fps, dtype=np.uint8)
## Using a flat index
d = target_products_fps_np.shape[1]                           # dimension
# original_dim = target_products_fps_np.shape[1]

index_flat = faiss.IndexBinaryFlat(d*8)  # build a flat (CPU) index
index_flat.add(target_products_fps_np)
k = 10
target_query_index = 3
query = np.array(target_products_fps_np[target_query_index])
D, I = index_flat.search(np.expand_dims(query, axis=0), k+1)
print('D:', D)
print('I:', I)
topk_indices = I[0]
topk_distances = D[0]

from rdkit.Chem import Draw

def draw_molecule_to_png(smi, file_name):
    '''
    use RDKit to draw molecule & save as [file_name](should be png)
    input:  1. smi: str, SMILES
            2. file_name : str, target file path
    '''
    mol = Chem.MolFromSmiles(smi)
    Draw.MolToFile(mol, file_name, size=(350, 300))

original_mol = target_products[target_query_index]
draw_molecule_to_png(original_mol, './test_find_knear2/original.png')
for i in range(k):
    index = topk_indices[i+1]
    mol = target_products[index]
    draw_molecule_to_png(mol, './test_find_knear2/rank%d.png' % (i+1))