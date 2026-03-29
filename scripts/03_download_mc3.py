import requests
import os
import time

DATA_DIR = "data/mafs"
os.makedirs(DATA_DIR, exist_ok=True)

TCGA_MAFS = {
    "TCGA-ACC":  "71cf7f4c-5b56-4425-b797-82bc9419a5d7",
    "TCGA-BLCA": "3eba6e6a-edc7-4daf-9a96-f4a9bbc8aefe",
    "TCGA-BRCA": "205f7120-573f-4dc4-94e2-1fb25302ff30",
    "TCGA-CESC": "f1ffcea6-9fcd-43bf-a0b3-9dc74c924c58",
    "TCGA-CHOL": "46e92822-8430-437a-a670-6432e2eed044",
    "TCGA-COAD": "35fa9eee-3082-4cf7-8b15-1dafafafe31e",
    "TCGA-ESCA": "ceeabfad-6e08-4bac-b9e5-14ce3fc4795f",
    "TCGA-GBM":  "99c228a0-02df-46f2-a095-68e72e53b159",
    "TCGA-HNSC": "ff591403-2ffa-4fa8-b133-8e3113ce8fa4",
    "TCGA-KICH": "79fbf0f9-9aac-4213-a8e2-8570b037e1ef",
    "TCGA-LGG":  "532d5274-fc66-486e-88f2-2ec6b147bb6c",
    "TCGA-LIHC": "226968c0-c961-4057-8119-0cd67ff5e4cb",
    "TCGA-LUAD": "d7718619-286f-47ec-9b8c-430fb1a01383",
    "TCGA-LUSC": "bf2ada35-4f78-44f0-9369-2a6f1e5824e1",
    "TCGA-OV":   "f424ca61-d189-4eb1-bf0d-02e2f44f1ea6",
    "TCGA-PAAD": "8c2214bf-daa9-4ea5-bf48-8c0cb0eba599",
    "TCGA-PRAD": "9ed651fe-63d8-4059-94a3-2a16639877f3",
    "TCGA-READ": "b77415e6-52ca-47ca-aa01-5ddb92f94e0f",
    "TCGA-SARC": "d66a02d0-5e7d-4a98-999c-4285514543ad",
    "TCGA-SKCM": "75047559-01f2-4146-95a8-ae679491df62",
    "TCGA-STAD": "5209c916-4975-4b04-ae57-457045098b4e",
    "TCGA-UCEC": "411a606d-f354-4561-b42c-eb1bc155c87d",
    "TCGA-UCS":  "48a4e0a8-d3cc-4b95-bd23-014a50c757e3",
}

BASE_URL = "https://api.gdc.cancer.gov/data/"

def download_file(file_id, output_path, max_retries=3):
    for attempt in range(max_retries):
        try:
            r = requests.get(BASE_URL + file_id, stream=True, timeout=60)
            if r.status_code == 200:
                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
                return True
            else:
                return False
        except Exception:
            if attempt < max_retries - 1:
                if os.path.exists(output_path):
                    os.remove(output_path)
                time.sleep(5)
            else:
                return False

print(f"Downloading {len(TCGA_MAFS)} cancer type MAFs...\n")
failed = []

for cancer_type, file_id in TCGA_MAFS.items():
    output_path = os.path.join(DATA_DIR, f"{cancer_type}.maf.gz")
    if os.path.exists(output_path):
        print(f"  [SKIP] {cancer_type}")
        continue
    print(f"  Downloading {cancer_type}...", end=" ", flush=True)
    success = download_file(file_id, output_path)
    if success:
        size = os.path.getsize(output_path) / 1e6
        print(f"✓ {size:.1f} MB")
    else:
        print("✗ Failed")
        failed.append(cancer_type)

print(f"\nDone! Failed: {failed if failed else 'none'}")