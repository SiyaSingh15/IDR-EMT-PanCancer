import requests
import os
import time

DATA_DIR = "data/expression"
os.makedirs(DATA_DIR, exist_ok=True)

# Direct download URL for TCGA PanCancer Atlas expression matrix
# Batch-effect corrected, log2(RSEM+1), all ~10,500 TCGA primary tumours
# Source: UCSC Xena pancanatlas hub
URL = "https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
OUT = "data/expression/tcga_pancan_expression.gz"

if os.path.exists(OUT):
    print(f"Already exists: {OUT} ({os.path.getsize(OUT)/1e6:.0f} MB)")
else:
    print("Downloading TCGA PanCancer expression matrix...")
    print(f"URL: {URL}\n")

    for attempt in range(3):
        try:
            r = requests.get(URL, stream=True, timeout=180)
            print(f"Status: {r.status_code}")
            if r.status_code == 200:
                total = 0
                with open(OUT, "wb") as f:
                    for chunk in r.iter_content(chunk_size=4*1024*1024):
                        f.write(chunk)
                        total += len(chunk)
                        print(f"  {total/1e6:.0f} MB", end="\r")
                print(f"\nDone: {os.path.getsize(OUT)/1e6:.0f} MB")
                break
            else:
                print(f"Failed: HTTP {r.status_code}")
                print(f"Headers: {dict(r.headers)}")
        except Exception as e:
            print(f"Attempt {attempt+1} failed: {e}")
            if attempt < 2:
                time.sleep(5)