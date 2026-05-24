import os
import zipfile
import urllib.request

EOS_URLS = {
    "APR": "https://compose.obspm.fr/download//1D/NS/Classical/APR/eos.zip",
    "SLY4": "https://compose.obspm.fr/download//1D/NS/Skyrme/SLY4/eos.zip",
    "DD2": "https://compose.obspm.fr/download//1D/NS/RMF/DD2/eos.zip",
}

BASE_DIR = "data"
os.makedirs(BASE_DIR, exist_ok=True)

for eos_name, url in EOS_URLS.items():
    out_dir = os.path.join(BASE_DIR, eos_name)
    zip_path = os.path.join(BASE_DIR, f"{eos_name}.zip")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 70)
    print(f"Downloading {eos_name}")
    print(f"Source: {url}")

    urllib.request.urlretrieve(url, zip_path)

    print(f"Extracting {zip_path}")
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(out_dir)

    print(f"{eos_name} files:")
    for f in sorted(os.listdir(out_dir)):
        print("  ", f)

print("=" * 70)
print("All EOS datasets downloaded and extracted successfully.")
