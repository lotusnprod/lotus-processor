import pubchempy as pcp

c = pcp.Compound.from_cid(5090)


print(c.molecular_formula)
print(c.molecular_weight)
print(c.isomeric_smiles)
print(c.xlogp)
print(c.iupac_name)
print(c.synonyms)
