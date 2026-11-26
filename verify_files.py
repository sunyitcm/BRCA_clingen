"""Check the exact release members and their contents; unrelated user files are ignored."""
from hashlib import sha256
from pathlib import Path
import re

MANIFEST = "MANIFEST.sha256"


def release_bytes(root):
    root = Path(root).resolve()
    mf = root / MANIFEST
    if mf.is_symlink() or not mf.is_file():
        raise ValueError("找不到发布清单")
    content = {MANIFEST: mf.read_bytes()}
    for line in content[MANIFEST].decode("utf-8").splitlines():
        m = re.fullmatch(r"([a-f0-9]{64})  ([A-Za-z0-9_./-]+)", line)
        if not m:
            raise ValueError("发布清单格式错误")
        digest, name = m.groups()
        parts = name.split("/")
        if any(s in ("", ".", "..", ".git") for s in parts) or name in content:
            raise ValueError("发布清单路径冲突：" + name)
        path = root
        for part in parts:
            path = path / part
            if path.is_symlink():
                raise ValueError("发布文件不能是链接：" + name)
        if not path.is_file():
            raise ValueError("发布文件缺失：" + name)
        data = path.read_bytes()
        if sha256(data).hexdigest() != digest:
            raise ValueError("文件校验失败：" + name)
        content[name] = data
    if len(content) == 1:
        raise ValueError("空发布清单")
    return content


def verify(root):
    return set(release_bytes(root)) - {MANIFEST}


if __name__ == "__main__":
    print("PASS · 清单内", len(verify(Path(__file__).resolve().parent)), "份文件逐一匹配；不检查额外的本地文件。")
