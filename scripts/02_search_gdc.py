import requests

search_url = "https://api.gdc.cancer.gov/files"

params = {
    "filters": '{"op":"and","content":[{"op":"=","content":{"field":"data_format","value":"MAF"}},{"op":"=","content":{"field":"access","value":"open"}},{"op":"=","content":{"field":"data_type","value":"Masked Somatic Mutation"}},{"op":"in","content":{"field":"cases.project.program.name","value":["TCGA"]}}]}',
    "fields": "file_id,file_name,file_size,cases.project.project_id",
    "size": "500",
    "sort": "file_size:desc",
    "expand": "cases.project"
}

r = requests.get(search_url, params=params).json()
hits = r["data"]["hits"]

print(f"Total files found: {len(hits)}\n")

# Getting unique projects and one representative file per project
projects = {}
for h in hits:
    try:
        proj = h["cases"][0]["project"]["project_id"]
        if proj not in projects:
            projects[proj] = {
                "file_id": h["file_id"],
                "file_name": h["file_name"],
                "size_mb": round(h["file_size"] / 1e6, 1)
            }
    except:
        pass

print(f"Unique TCGA cancer types: {len(projects)}\n")
for proj, info in sorted(projects.items()):
    print(f"  {proj}: {info['file_id']} | {info['size_mb']} MB")
