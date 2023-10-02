# This is horribly dirty…
# Also it is far from capturing everything


import glob
import pandas as pd
import re
import zipfile
from lxml import etree

problem_files = []
data = {}


def outstatus(message, color='red'):
    """Used to color the output for the terminal"""
    OK = '\033[92m'
    # WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    if color == 'red':
        return ('{} {} {}'.format(FAIL, message, ENDC))
    return ('{} {} {}'.format(OK, message, ENDC))


def treat_file(content):
    """We use the XML parser to read the string 'content' containing the xml
    file content"""
    try:
        tree = etree.XML(content)
        # We call the subsuck function that will extract data (not everything) from
        # the XML
        inchi = subsuck(data, tree)
        # If the inchi is not here, add it to a list of problematic files (we don't
        # do anything with that yet)
        if inchi is False:
            problem_files.append(filename)
        return inchi
    except:
        print("error")


def remove_tags(raw_html):
    """Remove some html from the xml string. This should be done differently, they
    add a lot of meaning with bold, and stuff like that."""
    # This function sucks really hard, it probably destroy a lot of things in
    # this file

    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext


def subsuck(data, element, parname=''):
    """The monster that suck the XML file and adds to the global variable data,
    this is horribly dirty, we should use classes and so on. But using XSLT
    would also make things easier, so I'm not rewriting that now.
    """
    inchi = False
    # Grab the compound code
    cid = element.attrib['chapman-hall-number']
    # If it is not in the data, create room for it
    if cid not in data:
        data[cid] = {}
    # Store it also as a column (avoids index columns)
    data[cid]['CRC_Number'] = cid
    parname += remove_tags(element.attrib['name'])
    data[cid]['Molecule_Name'] = parname
    # get is useful here to return an empty string if there is no element with
    # that name
    data[cid]['Molecule_Formula'] = remove_tags(element.attrib.get(
        'molecular-formula',
        ''))
    data[cid]['Molecule_Weight'] = element.attrib.get(
        'min_molecular-weight',
        '')
    data[cid]['Accurate_Mass'] = element.attrib.get('accurate-mass', '')
    data[cid]['CAS_Number'] = element.attrib.get('cas-number', '')
    data[cid]['Biological_Source'] = ''
    data[cid]['Compound_Types'] = ''
    data[cid]['Optical_Rotation'] = ''
    data[cid]['Biological_Use'] = ''
    data[cid]['Toxicity'] = ''
    data[cid]['Use'] = ''
    # Here we grab all the child elements that describe different parts of this
    # compound
    for child in element:
        if child.tag == 'compound-type':
            if data[cid]['Compound_Types'] != ' ':
                data[cid]['Compound_Types'] += ' '
            data[cid]['Compound_Types'] += remove_tags(child.attrib['type'])
        if child.tag == 'biosoc-sp2000':
            if data[cid]['Biological_Source'] != '':
                data[cid]['Biological_Source'] += ' '
            data[cid]['Biological_Source'] += remove_tags(child.attrib['text'])
        if child.tag == 'optical-rotation':
            if data[cid]['Optical_Rotation'] != '':
                data[cid]['Optical_Rotation'] += ' '
            data[cid]['Optical_Rotation'] += child.attrib['optical-rotation']
        if child.tag == 'biouse':
            if data[cid]['Biological_Use'] != '':
                data[cid]['Biological_Use'] += ' '
            data[cid]['Biological_Use'] += remove_tags(child.attrib['text'])
        if child.tag == 'hazard':
            if data[cid]['Toxicity'] != '':
                data[cid]['Toxicity'] += ' '
            data[cid]['Toxicity'] += remove_tags(child.attrib['text'])
        if child.tag == 'use':
            if data[cid]['Use'] != '':
                data[cid]['Use'] += ' '
            data[cid]['Use'] += remove_tags(child.attrib['text'])
        if child.tag == 'inchi':
            inchi = True
            chi = remove_tags(child.attrib['inchi'])
            # If the Inchi is not broken, lets grab it and the inchikey with it
            if (
                    'Unrecognized' not in chi and
                    'error' not in chi and
                    'Uninterpretable' not in chi):
                data[cid]['MolfileName'] = chi
                data[cid]['InChIKey'] = remove_tags(child.attrib['inchi-key'])
        # If there are variants of this molecule, also run that on them note
        # that we are passing parname, that will get added before the molecule
        # name because their "variants" just give the suffix (like glycoside)
        if child.tag == 'variant' or child.tag == 'derivative':
            inchi |= subsuck(data, child, parname)
    return inchi


# The main part of the program


print('DNP Cleaner')
print('- Treating jar files')
# We grab all the jar files
for filename in glob.glob('../data/external/dbSource/dnp/31_2/xmlfiles/*.jar'):
    count = 0
    # We open the zip files
    with zipfile.ZipFile(filename, 'r') as zip:
        status = True
        success = 0
        fail = 0
        num_files = len(zip.filelist)
        # For each file in the zip
        for file in zip.filelist:
            # If the file is a cache file we do something with it
            if '.cache' in file.filename:
                output = treat_file(zip.open(file.filename).read())
                if output is True:
                    status = status & output
                else:
                    status = status
                if output is True:
                    success += 1
                else:
                    # We count as failed any file that has no inchi
                    # they are not "really" failed
                    fail += 1
            # A dirty trick to put progress bars
            if (
                    int(100.0 * (count) / num_files) % 10 == 0 or
                    count >= num_files - 1):
                print('\r - {}  {}'.format(
                    filename.split("/")[2],
                    '#' * int(100.0 * (count + 1) / num_files)), end='')
            count += 1
        if status is True:
            print(outstatus(' OK for {} files'.format(success), 'green'))
        else:
            print(outstatus(' NOT OK for {}/{} files '.format(fail,
                                                              success), 'red'))

print("- Generating the CSV with the \"full set\" \
       (not everything is extracted yet)")
# We generate a panda dataframe with the table and make a csv out of it
df = pd.DataFrame.from_dict(list(data.values()))
order_rows = ["InChIKey", "MolfileName", "CRC_Number", "Molecule_Name",
              "Molecule_Formula", "Molecule_Weight", "Accurate_Mass",
              "CAS_Number", "Compound_Types", "Biological_Source",
              "Optical_Rotation", "Biological_Use", "Toxicity", "Use"]
# Any compound with an inchi gets exported in full_set
df[df.MolfileName.notnull()][order_rows].to_csv(
    "../data/external/dbSource/dnp/31_2/full_set.csv", index=False)

print("- Generating the CSV with the InChi")
order_rows_small = ["CRC_Number", "Molecule_Name",
                    "MolfileName"]
# This is just a file with less columns, not really useful
df[df.MolfileName.notnull()][order_rows_small].to_csv("../data/external/dbSource/dnp/31_2/id_inchi.csv",
                                                      index=False)
# If there is last year data
try:
    print("- Found last year data, generating the differences")
    # Load it
    olddf = pd.read_csv('last_year_data.csv')
    # For all inchikeys that were not here last year, add that to
    # new_structure.csv.
    # That means that any structure that changed will appear here
    # We should maybe do something more advanced here that will
    # tell what changed, if compounds get joined etc…
    df[~df.MolfileName.isin(olddf.MolfileName)][order_rows].to_csv(
        "new_structures.csv",
        index=False)
    # Get all the new CRC codes (that's the real new compounds)
    df[~df.CRC_Number.isin(olddf.CRC_Number)][order_rows].to_csv(
        "new_crc_codes.csv",
        index=False)
    # Get all the removed compounds
    # As they don't appear anymore in the files, we don't know why
    # they disappear
    olddf[~olddf.CRC_Number.isin(df.CRC_Number)][order_rows].to_csv(
        "removed_this_year.csv",
        index=False)
except:
    print("Weird error")
    pass
