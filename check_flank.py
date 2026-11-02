import sqlite3
con = sqlite3.connect("instance/ptm.db")
fa = open("verify_cache/P04637.fasta").read()
seq = "".join(fa.splitlines()[1:])
q = "select site, flanking, validation_status from ptm_sites where uniprot='P04637' and site in ('S6','S392')"
for site, fl, val in con.execute(q):
    pos = int(site[1:])
    print(site, "| record flank:", repr(fl), "| status:", val)
    print("    UniProt window (7 before, 7 after):", seq[max(0, pos - 8):pos + 7], "| residue", pos, "=", seq[pos - 1])