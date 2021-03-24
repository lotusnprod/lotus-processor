import pickle
import os
import numpy as np
import tmap as tm
import pandas as pd
import scipy.stats as ss
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from collections import Counter
import matplotlib.colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from time import sleep


def Top_N_classes(N=0, input_data=[], input_labels=[]):
    """Keep only the top N classes for the input classes and replace the following by 'Other', default N = 0 return input"""
    if N == 0:
        return input_data, input_labels
    else:
        top_input = [i for i, _ in Counter(input_data).most_common(N)]
        output_labels = [(15, "Other")]
        input_map = [15] * len(input_data)
        value = 1
        for i, name in input_labels:     
            if i in top_input:
                v = value
                if v == 15:
                    v = 0
                output_labels.append((v, name))
                input_map[i] = v
                value += 1
        output_data = [input_map[val] for _, val in enumerate(input_data)]
        return output_data, output_labels

def keep_only_given_class(to_keep = [], input_data=[], input_labels=[]):
    """Keep only the class specified in to_keep argument for mapping"""
    output_labels = []
    output_data = []
    i_to_keep = []
    for i, name in input_labels:
        if name in to_keep:
            output_labels.append(tuple([i, name]))
            i_to_keep.append(i)
        else:
            output_labels.append(tuple([1, 'Other']))
    output_labels = list(set(output_labels))

    for val in input_data:
        if val in i_to_keep:
            output_data.append(val)
        else:
            output_data.append(1)
    return output_data, output_labels

########################################################
### PART 1: LOAD Lotus DB and transform for plotting ###
########################################################

# Load lotus_db
df_meta = pd.read_csv('210223_frozen_metadata.csv.gz', sep=",")

# Fill NaN with 'Unknown'
values = {'organism_taxonomy_02kingdom': 'Not attributed (Bacteria and Algae)', 'organism_taxonomy_03phylum': 'Unknown', 'organism_taxonomy_04class':'Unknown',
'organism_taxonomy_05order': 'Unknown', 'organism_taxonomy_06family': 'Unknown', 'organism_taxonomy_08genus': 'Unknown','organism_taxonomy_09species': 'Unknown'}
df_meta.fillna(value=values, inplace=True)

# Keep only 1 occurence of structure - biosource pair
df_meta = df_meta.drop_duplicates(subset=['structure_wikidata', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum','organism_taxonomy_04class',
'organism_taxonomy_05order','organism_taxonomy_06family','organism_taxonomy_08genus', 'organism_taxonomy_09species'])

# Group by structure and return for each taxonomical level the most frequent taxa, its occurence and the number of unique taxa
df_gb = df_meta.groupby('structure_wikidata').agg(
    {   'structure_inchi':'first',
        'structure_inchikey': 'first',
        'structure_smiles': 'first',
        'structure_taxonomy_npclassifier_01pathway': 'first',
        'structure_taxonomy_npclassifier_02superclass': 'first',
        'structure_taxonomy_npclassifier_03class': 'first',
        'structure_taxonomy_classyfire_01kingdom':'first',
        'structure_taxonomy_classyfire_02superclass':'first',
        'structure_taxonomy_classyfire_03class':'first',
        'structure_taxonomy_classyfire_04directparent':'first',
        'organism_taxonomy_02kingdom': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'count', 'nunique'],
        'organism_taxonomy_03phylum': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        'organism_taxonomy_04class': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        'organism_taxonomy_05order': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        'organism_taxonomy_06family': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        'organism_taxonomy_08genus': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        'organism_taxonomy_09species': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'nunique'],
        })

del(df_meta)

df_gb.columns = ['_'.join(col).strip() for col in df_gb.columns.values]

df_gb.rename(columns={"organism_taxonomy_02kingdom_count": "biosource_count"}, inplace=True)

taxo_cols = ['organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum','organism_taxonomy_04class',
'organism_taxonomy_05order','organism_taxonomy_06family','organism_taxonomy_08genus', 'organism_taxonomy_09species']

bins = [0, 2, 6, 21, np.inf]
names = ['1 unique biosource', '2-5 unique biosources', '6-20 unique biosources', '>20 unique biosources']

for col in taxo_cols:
    df_gb[col+'_<lambda_0>'] = df_gb[col+'_<lambda_0>'].astype(str)
    df_gb[col+'_specificity'] = df_gb[col+'_<lambda_1>'] / df_gb['biosource_count']
    df_gb.drop([col+'_<lambda_1>'], axis=1, inplace=True)
    df_gb[col+'_nunique_cat'] = pd.cut(df_gb[col+'_nunique'], bins, labels=names, right = False)

names = ['1 biosource', '2-5 biosources', '6-20 biosources', '>20 biosources']
df_gb['biosource_count_cat'] = pd.cut(df_gb.biosource_count, bins, labels=names, right = False)

df_gb.reset_index(inplace=True)

df_gb["id_wikidata"] = df_gb['structure_wikidata'].str[31:]

values = {'structure_taxonomy_npclassifier_01pathway_first': 'Unknown',
    'structure_taxonomy_npclassifier_02superclass_first': 'Unknown',
    'structure_taxonomy_npclassifier_03class_first': 'Unknown',
    'structure_taxonomy_classyfire_01kingdom_first': 'Unknown',
    'structure_taxonomy_classyfire_02superclass_first': 'Unknown',
    'structure_taxonomy_classyfire_03class_first':'Unknown',
    'structure_taxonomy_classyfire_04directparent_first': 'Unknown'}
df_gb.fillna(value=values, inplace=True)

# Generating mean specificity scoreb by taxonomical level and class
dic_mean_specificity = {
'organism_taxonomy_02kingdom_specificity': {'min' : 0.9},
'organism_taxonomy_06family_specificity': {'min' : 0.8},
'organism_taxonomy_08genus_specificity': {'min' : 0.8},
}
structure_taxo = ['structure_taxonomy_npclassifier_01pathway_first', 'structure_taxonomy_npclassifier_02superclass_first', 'structure_taxonomy_npclassifier_03class_first']

for dic in dic_mean_specificity:
    for level in structure_taxo:
        df_gb['mean'] = df_gb.groupby([level])[str(dic)].transform('mean')
        df_gb.loc[df_gb[level] == 'Unknown', 'mean'] = df_gb["mean"].mean()
        df_gb.loc[df_gb['mean'] < dic_mean_specificity[dic]['min'], 'mean'] = dic_mean_specificity[dic]['min']
        dic_mean_specificity[dic][level] = df_gb['mean']
        df_gb.drop('mean', axis=1, inplace=True)

# Generating class for plotting
dic_categories = {
'structure_taxonomy_npclassifier_01pathway_first': {'Ncat' : 19}, 
"structure_taxonomy_npclassifier_02superclass_first" : {'Ncat' : 19},
"structure_taxonomy_npclassifier_03class_first": { 'Ncat' : 19},
"structure_taxonomy_classyfire_01kingdom_first": { 'Ncat' : 0},
"structure_taxonomy_classyfire_02superclass_first": { 'Ncat' : 19},
"structure_taxonomy_classyfire_03class_first": {'Ncat' : 19},
"structure_taxonomy_classyfire_04directparent_first": { 'Ncat' : 19},
"organism_taxonomy_02kingdom_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_03phylum_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_04class_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_05order_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_06family_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_08genus_<lambda_0>": { 'Ncat' : 19},
"organism_taxonomy_09species_<lambda_0>": { 'Ncat' :19},
"biosource_count_cat": { 'Ncat' : 0},
"organism_taxonomy_02kingdom_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_03phylum_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_04class_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_05order_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_06family_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_08genus_nunique_cat": { 'Ncat' : 0},
"organism_taxonomy_09species_nunique_cat": { 'Ncat' : 0}
}

for dic in dic_categories:
    labels, data = Faerun.create_categories(df_gb[str(dic)])
    dic_categories[dic]['data'], dic_categories[dic]['labels'] = Top_N_classes(dic_categories[dic]['Ncat'], data, labels)

# Create a categories to map only a given taxononomical family

labels, data = Faerun.create_categories(df_gb['organism_taxonomy_06family_<lambda_0>'])
simaroubaceae_data, simaroubaceae_labels = keep_only_given_class(['Simaroubaceae'], data, labels)

# Generating colormaps for plotting
cmap = mcolors.ListedColormap(["gainsboro", "peachpuff", "salmon", "tomato"])
cmap2 = mcolors.ListedColormap(["gainsboro", "tomato"])
cmap3 = mcolors.ListedColormap(["lightgray", "lightgray","lightgray", "lightgray", "lightgray","lightgray", "lightgray", "firebrick"])

# Generate a label column
df_gb["labels"] = (
    df_gb["structure_smiles_first"]
    + '__<a target="_blank" href="'
    + df_gb["structure_wikidata"]
    + '">'
    + df_gb["id_wikidata"]
    + '__'
    + df_gb['structure_inchikey_first']
    + '__'
    + df_gb["structure_taxonomy_npclassifier_03class_first"]
    + "</a>"
)


########################################################
# OPTION 1: START FROM SCRATCH AND GENERATE LSH FOREST AND TMAP FROM SMILES, PLUS CHEMICAL DESCRIPTORS
########################################################

enc = MHFPEncoder(1024)
lf = tm.LSHForest(1024, 64)

fps = []
hac = []
c_frac = []
ring_atom_frac = []
largest_ring_size = []

n_iter = len(df_gb)
for i, row in df_gb.iterrows():
    j = (i + 1) / n_iter
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(50 * j):{50}s}] {round((100 * j),2)}%")
    sys.stdout.flush()
    sleep(0.05)

    mol = AllChem.MolFromInchi(row["structure_inchi_first"])
    atoms = mol.GetAtoms()
    size = mol.GetNumHeavyAtoms()
    n_c = 0
    n_ring_atoms = 0
    for atom in atoms:
        if atom.IsInRing():
            n_ring_atoms += 1
        if atom.GetSymbol().lower() == "c":
            n_c += 1

    c_frac.append(n_c / size if size != 0 else 0)

    ring_atom_frac.append(n_ring_atoms / size if size != 0 else 0)
        
    sssr = AllChem.GetSymmSSSR(mol)
    if len(sssr) > 0:
        largest_ring_size.append(max([len(s) for s in sssr]))
    else:
        largest_ring_size.append(0)
    hac.append(size)
    fps.append(tm.VectorUint(enc.encode_mol(mol)))

lf.batch_add(fps)
lf.index()

#Store lsh forest and structure metadata
lf.store("210312_lotus_db.dat")
with open("210312_lotus_db.pickle", "wb+") as f:
    pickle.dump(
        (hac, c_frac, ring_atom_frac, largest_ring_size),
        f,
        protocol=pickle.HIGHEST_PROTOCOL,
    )
    
# tmap configuration
cfg = tm.LayoutConfiguration()
cfg.node_size = 1 / 50
cfg.mmm_repeats = 2
cfg.sl_repeats = 2

# tmap generation
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
del(lf)


########################################################
# OPTION 2: LOAD PRE_COMPUTED LSH FOREST AND CHEMICAL DESCRIPTORS
########################################################

lf = tm.LSHForest(1024, 64)
lf.restore("210312_lotus_db.dat") # Version "210312_lotus_db.dat" contains 270'336 resulting from 210223_frozen_metadata groupby(structure_wikidata)

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
        open("210312_lotus_db.pickle", "rb")
)
c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

# tmap configuration
cfg = tm.LayoutConfiguration()
cfg.node_size = 1 / 50
cfg.mmm_repeats = 2
cfg.sl_repeats = 2

# tmap generation
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg) 
del(lf)

########################################################
# OPTION 3: LOAD PRE-COMPUTED COORDINATES, SOURCES AND TARGETS AND CHEMICAL DESCRIPTORS
########################################################

x, y, s, t = pickle.load(open("coords_210312_lotus_db.dat", "rb")) # Version "coords_210312_lotus_db.dat" contains 270'336 resulting from 210223_frozen_metadata groupby(structure_wikidata)

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
        open("210312_lotus_db.pickle", "rb")
)
c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

########################################################
# COMMON PART: PLOT THE TMAP
########################################################

# Plotting function
f = Faerun(view="front", coords=False, clear_color='#ffffff')
f.add_scatter(
    "lotus_db",
    {
        "x": x,
        "y": y,
        "c": [
            dic_categories['structure_taxonomy_npclassifier_01pathway_first']['data'],
            dic_categories['structure_taxonomy_npclassifier_02superclass_first']['data'],
            dic_categories['structure_taxonomy_npclassifier_03class_first']['data'],
            dic_categories['structure_taxonomy_classyfire_01kingdom_first']['data'],
            dic_categories['structure_taxonomy_classyfire_02superclass_first']['data'],
            dic_categories['structure_taxonomy_classyfire_03class_first']['data'],
            dic_categories['structure_taxonomy_classyfire_04directparent_first']['data'],
            dic_categories['organism_taxonomy_02kingdom_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_03phylum_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_04class_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_05order_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_06family_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_08genus_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_09species_<lambda_0>']['data'],
            dic_categories['biosource_count_cat']['data'],
            dic_categories['organism_taxonomy_02kingdom_nunique_cat']['data'],
            dic_categories['organism_taxonomy_03phylum_nunique_cat']['data'],
            dic_categories['organism_taxonomy_04class_nunique_cat']['data'],
            dic_categories['organism_taxonomy_05order_nunique_cat']['data'],
            dic_categories['organism_taxonomy_06family_nunique_cat']['data'],
            dic_categories['organism_taxonomy_08genus_nunique_cat']['data'],
            dic_categories['organism_taxonomy_09species_nunique_cat']['data'],
            simaroubaceae_data,
            dic_mean_specificity['organism_taxonomy_02kingdom_specificity']['structure_taxonomy_npclassifier_01pathway_first'],
            dic_mean_specificity['organism_taxonomy_02kingdom_specificity']['structure_taxonomy_npclassifier_02superclass_first'],
            dic_mean_specificity['organism_taxonomy_02kingdom_specificity']['structure_taxonomy_npclassifier_03class_first'],
            dic_mean_specificity['organism_taxonomy_06family_specificity']['structure_taxonomy_npclassifier_01pathway_first'],
            dic_mean_specificity['organism_taxonomy_06family_specificity']['structure_taxonomy_npclassifier_02superclass_first'],
            dic_mean_specificity['organism_taxonomy_06family_specificity']['structure_taxonomy_npclassifier_03class_first'],
            dic_mean_specificity['organism_taxonomy_08genus_specificity']['structure_taxonomy_npclassifier_01pathway_first'],
            dic_mean_specificity['organism_taxonomy_08genus_specificity']['structure_taxonomy_npclassifier_02superclass_first'],
            dic_mean_specificity['organism_taxonomy_08genus_specificity']['structure_taxonomy_npclassifier_03class_first'],
            hac,
            c_frak_ranked,
            ring_atom_frac,
            largest_ring_size,
        ],
        "labels": df_gb["labels"],
     },
    shader="smoothCircle",
    point_scale=2.0,
    max_point_size=10,
    legend_labels=[
        dic_categories['structure_taxonomy_npclassifier_01pathway_first']['labels'],
        dic_categories['structure_taxonomy_npclassifier_02superclass_first']['labels'],
        dic_categories['structure_taxonomy_npclassifier_03class_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_01kingdom_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_02superclass_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_03class_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_04directparent_first']['labels'],
        dic_categories['organism_taxonomy_02kingdom_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_03phylum_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_04class_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_05order_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_06family_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_08genus_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_09species_<lambda_0>']['labels'],
        dic_categories['biosource_count_cat']['labels'],
        dic_categories['organism_taxonomy_02kingdom_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_03phylum_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_04class_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_05order_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_06family_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_08genus_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_09species_nunique_cat']['labels'],
        simaroubaceae_labels
    ],
    categorical=[True, True, True, True,  True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,True, True,
    False, False, False,False, False,False, False,False, False,False, False, False, False],
    colormap=["tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20",  cmap, cmap2, cmap, cmap, cmap, cmap, cmap, cmap, cmap3,
    "inferno","inferno","inferno","inferno","inferno","inferno","inferno","inferno","inferno","rainbow",  "rainbow", "rainbow","Blues"],
    series_title=[
        "chemo_np_classifier_pathway",
        "chemo_np_classifier_superclass",
        "chemo_np_classifier_class",
        "chemo_cf_kingdom",
        "chemo_cf_superclass",
        "chemo_cf_class",
        "chemo_cf_parent",
        "kingdom_bio",
        "phylum_bio",
        "class_bio",
        "order_bio",
        "family_bio",
        "genus_bio",
        "species_bio",
        "biosource_count",
        "kingdom_nunique",
        "phylum_nunique",
        "class_nunique",
        "order_nunique",
        "family_nunique",
        "genus_nunique",
        "species_nunique",
        "Family",
        "Pathway kingdom specificity mean",
        "Superclass kingdom specificity mean",
        "Class kingdom specificity mean",
        "Pathway family specificity mean",
        "Superclass family specificity mean",
        "Class family specificity mean",
        "Pathway genus specificity mean",
        "Superclass genus specificity mean",
        "Class genus specificity mean",
        "HAC",
        "C Frac",
        "Ring Atom Frac",
        "Largest Ring Size",
    ],
    has_legend=True,
)
f.add_tree("lotus_db_tree", {"from": s, "to": t}, point_helper="lotus_db", color = '#e6e6e6')
f.plot('lotus_db', template="smiles")