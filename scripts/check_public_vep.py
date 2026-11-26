# 保存为 check_public_vep.py；运行 python3 check_public_vep.py
import json
from pathlib import Path

for name, original in [('vep-vcf.json', '13 32380106 . G T . . .'),
                       ('vep-c.json', 'NM_000059.4:c.9217G>T')]:
    try:
        items = json.loads(Path(name).read_bytes())
    except (OSError, ValueError) as exc:
        raise SystemExit(f'{name}: 文件不可读或不是 JSON，请看 HTTP 状态和响应头') from exc
    if not (isinstance(items, list) and len(items) == 1 and isinstance(items[0], dict)):
        raise SystemExit(f'{name}: 未得到唯一输入结果')
    rec = items[0]
    if (rec.get('assembly_name'), rec.get('input')) != ('GRCh38', original):
        raise SystemExit(f'{name}: 组装或输入回显不符')
    hits = [t for t in (rec.get('transcript_consequences') or []) if isinstance(t, dict)
            and t.get('gene_symbol') == 'BRCA2' and t.get('transcript_id') == 'NM_000059.4'
            and t.get('hgvsc') == 'NM_000059.4:c.9217G>T']
    if len(hits) != 1:
        raise SystemExit(f'{name}: 没有唯一匹配的规范转录本')
    print(name, hits[0]['hgvsc'], hits[0].get('hgvsp') or '未返回蛋白描述')
