import requests
import os
import time

DATA_DIR = "data/mafs"
os.makedirs(DATA_DIR, exist_ok=True)

BASE_URL = "https://api.gdc.cancer.gov/"

# All 33 TCGA cancer types
TCGA_PROJECTS = [
    "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
    "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
    "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
    "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
    "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
    "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
    "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
]

def get_all_maf_ids_for_project(project_id):
    """Get ALL open-access MAF file IDs for a given TCGA project."""
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_format", "value": "MAF"}},
            {"op": "=", "content": {"field": "access", "value": "open"}},
            {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}}
        ]
    }

    import json
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size",
        "size": "500"
    }

    r = requests.get(BASE_URL + "files", params=params)
    hits = r.json()["data"]["hits"]
    return hits

def download_file(file_id, output_path, max_retries=3):
    for attempt in range(max_retries):
        try:
            r = requests.get(BASE_URL + f"data/{file_id}", stream=True, timeout=120)
            if r.status_code == 200:
                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
                return True
            else:
                print(f"    Status {r.status_code}")
                return False
        except Exception as e:
            if attempt < max_retries - 1:
                if os.path.exists(output_path):
                    os.remove(output_path)
                print(f"    Retry {attempt+2}/{max_retries} in 5s...")
                time.sleep(5)
            else:
                return False

print("Fetching all TCGA MAF files per cancer type...\n")

total_downloaded = 0
project_summary = {}

for project in TCGA_PROJECTS:
    print(f"\n{project}:")
    hits = get_all_maf_ids_for_project(project)

    if not hits:
        print(f"  No open-access MAF files found")
        project_summary[project] = 0
        continue

    print(f"  Found {len(hits)} MAF files")
    project_dir = os.path.join(DATA_DIR, project)
    os.makedirs(project_dir, exist_ok=True)

    downloaded = 0
    for hit in hits:
        file_id = hit["file_id"]
        fname = hit["file_name"]
        size_mb = hit["file_size"] / 1e6
        out_path = os.path.join(project_dir, fname)

        if os.path.exists(out_path):
            downloaded += 1
            continue

        print(f"  Downloading {fname} ({size_mb:.1f} MB)...", end=" ", flush=True)
        success = download_file(file_id, out_path)
        if success:
            print("✓")
            downloaded += 1
            total_downloaded += 1
        else:
            print("✗ Failed")

    project_summary[project] = downloaded
    print(f"  {downloaded}/{len(hits)} files ready")

print(f"\n=== SUMMARY ===")
for proj, count in project_summary.items():
    print(f"  {proj}: {count} files")
pr\int(f"\nTotal new downloads: {total_downloaded}")