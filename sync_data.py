
"""
Sync the real records from the original PTM-SiteDB search pages.

Run:
    python sync_data.py

It downloads the paginated search results from the original site,
extracts the table fields already displayed by PTM-SiteDB, and stores
them locally in instance/ptm.db. No records are invented.
"""
import re, sqlite3, time, datetime
from pathlib import Path
from urllib.parse import urlencode
import requests
from bs4 import BeautifulSoup

import sys
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE))
DB=BASE/"instance"/"ptm.db"
ROOT="http://93.127.214.36/ptmsitedb/pages/search.php"
HEADERS={"User-Agent":"PTM-SiteDB-Student-Project/1.0"}

def clean(s):
    return re.sub(r"\s+"," ",s or "").strip()

def parse_page(html, source_url):
    soup=BeautifulSoup(html,"html.parser")
    table=soup.find("table")
    if not table: return []
    out=[]
    for tr in table.select("tbody tr"):
        tds=tr.find_all("td")
        if len(tds)<10: continue
        gene=clean(tr.select_one(".td-gene").get_text(" ",strip=True) if tr.select_one(".td-gene") else tds[0].get_text(" ",strip=True))
        uni=clean(tds[1].get_text(" ",strip=True))
        ptm=clean(tds[2].get_text(" ",strip=True))
        site=clean(tr.select_one(".td-site").get_text(" ",strip=True) if tr.select_one(".td-site") else tds[3].get_text(" ",strip=True))
        flank=clean(tr.select_one(".td-flanking").get_text(" ",strip=True) if tr.select_one(".td-flanking") else tds[4].get_text(" ",strip=True))
        flank=re.sub(r"\s+","",flank)
        kinase=clean(tr.select_one(".td-kinase").get_text(" ",strip=True) if tr.select_one(".td-kinase") else tds[5].get_text(" ",strip=True))
        disease=clean(tds[6].get_text(" ",strip=True))
        evidence=clean(tds[7].get_text(" ",strip=True))
        confidence=clean(tds[8].get_text(" ",strip=True))
        pmid=clean(tds[9].get_text(" ",strip=True))
        m=re.search(r"(\d+)$",site)
        pos=int(m.group(1)) if m else None
        out.append((gene,uni,ptm,site,pos,flank,kinase,disease,evidence,confidence,pmid,source_url))
    return out

def sync():
    DB.parent.mkdir(exist_ok=True)
    conn=sqlite3.connect(DB)
    from services.dbutil import SCHEMA
    conn.executescript(SCHEMA)
    before=conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    total=0; pages=0; source_error=None
    for page in range(1,401):
        params={"page":page}
        # No ptm_type filter: preserve all PTM types available in the original site.
        url=ROOT+"?"+urlencode(params)
        try:
            r=requests.get(url,headers=HEADERS,timeout=20)
            r.raise_for_status()
        except requests.RequestException as exc:
            source_error=f"Could not connect to original PTM-SiteDB source on page {page}: {exc}"
            print("SOURCE CONNECTION ERROR:", source_error)
            break
        records=parse_page(r.text,url)
        if not records:
            break
        for rec in records:
            conn.execute("""INSERT OR IGNORE INTO ptm_sites
            (gene,uniprot,ptm_type,site,position,flanking,kinase,disease_context,evidence,confidence,pmid,source_url)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",rec)
        conn.commit()
        total += len(records); pages += 1
        print(f"page {page}: {len(records)} records")
        # The source shows up to 25 rows per page. A short page is the end.
        if len(records) < 25:
            break
        time.sleep(0.3)
    after=conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    conn.execute("INSERT INTO sync_log(synced_at,source,pages,records_parsed,records_new) VALUES (?,?,?,?,?)",
                 (datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None).isoformat(timespec="seconds") + "Z",ROOT,pages,total,after-before))
    conn.commit(); conn.close()
    if total==0:
        if source_error:
            print("WARNING: The source could not be reached. No substitute data were imported.")
        else:
            print("WARNING: 0 records parsed. The source page layout may have changed; the database was NOT filled with any substitute data.")
    print(f"Synced {total} displayed records into {DB}")
    return 2 if total == 0 else 0

if __name__=="__main__":
    raise SystemExit(sync())
