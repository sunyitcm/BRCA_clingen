"""Indexed local FASTA, pinned to an explicitly selected reference manifest."""
import hashlib
import json
import re
from pathlib import Path


def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Reference:
    def __init__(self, manifest, assembly):
        manifest = Path(manifest)
        raw = manifest.read_bytes()
        metadata = json.loads(raw)
        if metadata.get("assembly") != assembly:
            raise ValueError("FASTA 清单与输入组装不一致")
        self.assembly = assembly
        self.path = manifest.parent / metadata["fasta_path"]
        index = Path(str(self.path) + ".fai")
        for path, field in [(self.path, "fasta_SHA256"), (index, "fai_SHA256")]:
            if file_sha(path) != metadata.get(field):
                raise ValueError(f"{field} 与清单不符")
        self.provenance = {k: metadata.get(k) for k in
                           ("assembly", "fasta_SHA256", "fai_SHA256", "source_release")}
        self.provenance["manifest_sha256"] = hashlib.sha256(raw).hexdigest()
        self.entries = {}
        for line in index.read_text().splitlines():
            name, size, offset, bases, width = line.split("\t")[:5]
            if name not in ("13", "17", "chr13", "chr17"):
                continue
            chrom = name.removeprefix("chr")
            size, offset, bases, width = map(int, (size, offset, bases, width))
            if (chrom in self.entries or size != metadata.get("contigs", {}).get(name)
                    or offset < 1 or bases < 1 or width < bases or size < 1):
                raise ValueError("FAI 字典、行宽或染色体别名冲突")
            self.entries[chrom] = (size, offset, bases, width)

    def sequence(self, chrom, start, end):
        """Read [start,end) from the forward strand."""
        if chrom not in self.entries:
            raise ValueError("FASTA 缺少目标染色体")
        size, offset, bases, width = self.entries[chrom]
        if not 0 <= start <= end <= size:
            raise ValueError("参考区间越界")
        if start == end:
            return ""
        first = offset + start // bases * width + start % bases
        last = offset + (end-1) // bases * width + (end-1) % bases
        with self.path.open("rb") as handle:
            handle.seek(first)
            seq = b"".join(handle.read(last-first+1).splitlines()).upper()
        if len(seq) != end-start or not re.fullmatch(b"[ACGT]+", seq):
            raise ValueError("参考区间有未知碱基或索引不符")
        return seq.decode("ascii")

    def check(self, variant):
        if variant["assembly"] != self.assembly:
            raise ValueError("参考组装不一致")
        pos = variant["pos"]
        end = pos + len(variant["ref"]) - 1
        seq = self.sequence(variant["chrom"], pos-1, end)
        return {"id": f"chromosome:{self.assembly}:{variant['chrom']}:{pos}:{end}:1",
                "seq": seq, "source": "local_fasta", **self.provenance}

    def equivalent(self, chrom, left, right):
        """Verify both REF strings and compare alternate haplotypes in one window."""
        start = min(left[0], right[0])
        end = max(left[1], right[1])
        if end-start > 10000:
            return False
        genome = self.sequence(chrom, start, end)
        alternatives = []
        for a, b, ref, alt in (left, right):
            if b < a or b-a != len(ref) or genome[a-start:b-start] != ref:
                return False
            alternatives.append(genome[:a-start] + alt + genome[b-start:])
        return alternatives[0] == alternatives[1]
