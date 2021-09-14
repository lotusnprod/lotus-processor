# Antibase converter
# Converting from SDF to csv

# loading packages
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.rdmolfiles import SmilesWriter

# file path
my_sdf_file = '../data/external/dbSource/antimarin/antimarin0311_test.sdf'
my_smi_file = '../data/external/dbSource/antimarin/antimarinConverted.smi'
my_tsv_file = '../data/external/dbSource/antimarin/antimarinConverted.tsv.gz'

# converting
sdf_frame = PandasTools.LoadSDF(my_sdf_file,
                                smilesName=None,
                                embedProps=True,
                                molColName=None,
                                includeFingerprints=False)

# to get smiles
mols = [mol for mol in Chem.SDMolSupplier(my_sdf_file) if mol != None]

# make writer object with a file name.
writer = SmilesWriter(my_smi_file)
 
# SetProps method can set properties that will be written to files with SMILES.
writer.SetProps(['Mol_ID'])

# exporting
sdf_frame.to_csv(my_tsv_file,
                 compression='gzip', sep='\t')

# The way of writing molecules can perform common way.
for mol in mols:
    writer.write( mol )
writer.close()