# COCONUT converter
# Converting from SDF to csv

# loading packages
from rdkit.Chem import PandasTools

# file path
my_sdf_file = '../data/external/dbSource/coconut/COCONUT_DB.sdf'

# converting
sdf_frame = PandasTools.LoadSDF(my_sdf_file,
                                smilesName=None,
                                embedProps=True,
                                molColName=None,
                                includeFingerprints=False)

# exporting
sdf_frame.to_csv('../data/external/dbSource/coconut/coconutConverted.tsv.gz',
                 compression='gzip', sep='\t')
