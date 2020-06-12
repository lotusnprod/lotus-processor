# COCONUT converter
## Converting from SDF to csv

# loading packages
from rdkit import Chem
from rdkit.Chem import PandasTools
import zipfile

# file path
my_zip_file = zipfile.ZipFile(file = '../data/external/dbSource/COCONUT/0_initial_files/COCONUT.sdf.zip', mode = 'r')
my_sdf_file = my_zip_file.open('COCONUT.sdf')


# converting
sdf_frame = PandasTools.LoadSDF(my_sdf_file,
                                smilesName = None,
                                embedProps = True,
                                molColName = None,
                                includeFingerprints = False)

#exporting
sdf_frame.to_csv('../data/external/dbSource/COCONUT/0_initial_files/coconutCompiled.tsv.zip', compression = 'gzip', sep = '\t')
