import requests
import json
import gzip
from bs4 import BeautifulSoup
from multiprocessing.dummy import Pool  # This is a thread-based Pool
from multiprocessing import cpu_count

def tkm_process_row(row):
    return {
        "kind": row[0],
        "name": BeautifulSoup(row[1], "lxml").a.contents[0],
        "source": row[4],
        "relationship": row[5]
    }

def get_tkm_data(id_):
     page = requests.get(f"http://pharmdb-k.org/class/getRelate.php?q=2&category=TKM-Drug&id=10290825&type=drug&_search=false&rows=100000&page=1&sidx=Name_A&sord=asc")
     rows = json.loads(page.content)["rows"]
     return [tkm_process_row(row["cell"]) for row in rows]

def process_id(id_, retries=3):
    if (retries<0):
        print(f"FAILED: process_id {id_}")
        return {}
    page = requests.get(f"http://pharmdb-k.org/detail_view.php?category=Drug&nid={id_}")
    soup = BeautifulSoup(page.content, "lxml")
    table = soup.find(class_="table2")
    if table is None:
        return process_id(id_, retries-1)

    trs = table.find_all("tr")
    pubchem_cid = ""
    smiles = ""
    inchi = ""
    name = ""
    for tr in trs:
        if tr.th.contents[0] == "Compound Name":
            if (len(tr.td.contents)>0):
                name = tr.td.contents[0]
            else:
                return {}
        if tr.th.contents[0] == "PubChem CID":
            pubchem_cid = tr.td.a.contents[0]
        if tr.th.contents[0] == "Canonical SMILES":
            smiles = tr.td.contents[0].strip()
        if tr.th.contents[0] == "InChi":
            inchi = tr.td.contents[0].strip()

    tkm_data = get_tkm_data(pubchem_cid)
    return {"name": name, "pubchem_cid": pubchem_cid, "smiles": smiles, "inchi": inchi, "tkm_data": tkm_data}

page_number = 1
safety_max_page = 3000  # Increase this number if there are more pages

NUM_WORKERS = 10
pool = Pool(NUM_WORKERS)


while page_number < safety_max_page:
    out = []
    page = requests.get(f"http://pharmdb-k.org/pharmdb_search.php?category=All&sval=&GotoPage={page_number}")
    print(page_number)
    soup = BeautifulSoup(page.content, "lxml")
    ids = [block.find_all("td")[1].a["href"].split("=")[2] for block in soup.find(class_="table2").tbody.find_all("tr")]

    # Process in parallel and remove problematic entries
    out += [val for val in pool.imap(process_id, ids) if val != {}]

    # Check if we are on the last page
    last_element_in_page_buttons = soup.find(class_="tblPage").find_all("span")[-1]
    if last_element_in_page_buttons["class"][0] != "pageBt":
        print("On Last page")
        break

    page_number += 1
    with gzip.open(f"pharmdbk_page_{page_number}.json.gz", "wt") as f:
        f.write(json.dumps(out))

