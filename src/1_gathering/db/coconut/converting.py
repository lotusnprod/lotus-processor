# COCONUT converter
## Converting from SDF to csv

# loading packages
from rdkit import Chem
from rdkit.Chem import PandasTools
import zipfile

# file path
# my_zip_file = zipfile.ZipFile(file = '../data/external/dbSource/coconut/COCONUT.sdf.zip', mode = 'r')
# my_sdf_file = my_zip_file.open('COCONUT.sdf')

my_sdf_file = '../data/external/dbSource/coconut/COCONUT_DB.sdf'


# converting
sdf_frame = PandasTools.LoadSDF(my_sdf_file,
                                smilesName = None,
                                embedProps = True,
                                molColName = None,
                                includeFingerprints = False)

# sexporting
sdf_frame.to_csv('../data/external/dbSource/coconut/coconutConverted.tsv.zip', compression = 'gzip', sep = '\t')
