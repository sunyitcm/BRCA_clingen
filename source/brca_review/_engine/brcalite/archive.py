"""Read recorded HTTP evidence only when it matches the request exactly."""
import json
import shutil
from .net import getj, sha


def replay_fetch(previous, root, label, url, payload=None, kind="json"):
    source = previous / "http" / label
    receipt = getj(source / "receipt.json")
    body = None if payload is None else json.dumps(payload, separators=(",", ":")).encode()
    expected = (url, "GET" if body is None else "POST", None if body is None else sha(body))
    actual = tuple(receipt.get(k) for k in ("url", "method", "request_sha256"))
    if actual != expected or receipt.get("label") != label:
        raise ValueError(f"历史请求与此次查询不符 {label}")
    if body is not None and (source / "request.json").read_bytes() != body:
        raise ValueError(f"历史请求体不符 {label}")
    files = ["receipt.json"] + (["request.json"] if body is not None else [])
    attempts = receipt["attempts"]
    deferred = not attempts and receipt.get("state") == "host_cooldown" and receipt.get("retry_at_utc")
    if (not attempts and not deferred) or bool(attempts and attempts[-1].get("ok")) != receipt["ok"]:
        raise ValueError(f"历史请求状态不完整 {label}")
    result = None
    for i, attempt in enumerate(attempts, 1):
        name = attempt.get("body_file")
        if not name and not attempt.get("ok"):
            continue
        if name != ("body" if attempt.get("ok") else f"failed-{i}.body"):
            raise ValueError(f"历史响应文件不符 {label}")
        raw = (source / name).read_bytes()
        if sha(raw) != attempt["body_sha256"] or len(raw) != attempt["bytes"]:
            raise ValueError(f"历史响应校验失败 {label}")
        files.append(name)
        if attempt.get("ok"):
            result = raw
    target = root / label
    target.mkdir(parents=True)
    for name in files:
        shutil.copyfile(source / name, target / name)
    return result
