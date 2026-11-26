"""Build the RC asset and its checksum from the approved manifest, without network access."""
import argparse
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verify_files import release_bytes

STEM = "brca-review-desk-0.2.0-rc.8"
FIXED_ZIP_DATE = (1980, 1, 1, 0, 0, 0)  # Format epoch, unrelated to a release or capture.


def build(root, out):
    files = release_bytes(root)
    target = BytesIO()
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(files.items()):
            entry = zipfile.ZipInfo(f"{STEM}/{name}", FIXED_ZIP_DATE)
            entry.create_system = 3
            entry.external_attr = (0o100755 if name == "scripts/demo_vep_curl.sh" else 0o100644) << 16
            archive.writestr(entry, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    raw = target.getvalue()
    name = STEM + ".zip"
    digest = sha256(raw).hexdigest()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    outputs = ((out / name, raw), (out / "SHA256SUMS.txt", (digest + "  " + name + "\n").encode("ascii")))
    for dest, data in outputs:
        if dest.is_symlink() or (dest.exists() and (not dest.is_file() or dest.read_bytes() != data)):
            raise ValueError("输出目录已有不同的同名文件；请改用空目录。")
    for dest, data in outputs:
        if not dest.exists():
            # No filesystem owner, resource fork, time or local Git metadata enters the archive.
            with dest.open("xb") as f:
                f.write(data)
    return name, digest


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output-dir", type=Path, default=Path("dist"))
    args = cli.parse_args()
    name, digest = build(ROOT, args.output_dir)
    print(name, "SHA-256", digest)
