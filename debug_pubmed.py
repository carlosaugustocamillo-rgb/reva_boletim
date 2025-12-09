import os
from Bio import Entrez
import time
from dotenv import load_dotenv

load_dotenv()

ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL")
if not ENTREZ_EMAIL:
    print("❌ ENTREZ_EMAIL não configurado.")
    exit(1)

Entrez.email = ENTREZ_EMAIL

def debug_pubmed_search():
    query = '("cancer"[MeSH Terms] OR "oncology"[All Fields]) AND ("rehabilitation"[MeSH Terms] OR "physical therapy modalities"[MeSH Terms]) AND ("2025/11/26"[Date - Publication] : "2025/12/03"[Date - Publication])'
    
    print(f"🔎 Buscando por: {query}")
    start = time.time()
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        ids = record["IdList"]
        
        end = time.time()
        print(f"✅ Sucesso! Encontrados: {count} artigos.")
        print(f"⏱️ Tempo: {end - start:.2f}s")
        print(f"IDs: {ids}")
        
    except Exception as e:
        print(f"❌ Erro na busca: {e}")

if __name__ == "__main__":
    debug_pubmed_search()
