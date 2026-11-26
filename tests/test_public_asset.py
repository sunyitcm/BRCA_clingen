"""Release regressions: published examples, literal response hashes and deterministic archive."""
import base64
import hashlib
import importlib.util
import json
from html import unescape
from pathlib import Path
import re
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source"))
from brca_review.desk import check_original, render
from brca_review.examples import examples
from brca_review.model import compare

s = importlib.util.spec_from_file_location("build_public_asset", ROOT / "scripts/build_release.py")
builder = importlib.util.module_from_spec(s)
s.loader.exec_module(builder)


class PublicTeachingAndRelease(unittest.TestCase):
    def test_teaching_data_and_raw_reports_preserve_identity_without_capture_time(self):
        cases = examples()
        self.assertEqual(len(cases), 5)
        self.assertEqual(len({c['id'] for c in cases}), 5)
        for c in cases:
            self.assertTrue(c['id'].startswith('demo-public-rc8-'))
            self.assertIsNone(c['utc'])
            self.assertIsNone(c['report']['generated_utc'])
            self.assertFalse(c['raw'] or c['blobs'])
            check_original(c)
            raw = json.loads(base64.b64decode(c['raw_report'], validate=True))
            self.assertIsNone(raw['generated_utc'])
            self.assertEqual((c['caid'], c['assembly']), (raw['caid'], raw['input']['assembly']))
        c = next(c for c in cases if c['id'].endswith('brca2-38'))
        self.assertTrue(c['report']['erepo']['published'])
        self.assertIn('公开教学示例', compare(c, c)['message'])

    def test_bundled_html_is_rebuilt_from_the_public_inputs(self):
        expected = render(examples(), selected='demo-public-rc8-brca2-38')
        for path in (ROOT/'docs/snapshot.html', ROOT/'source/brca_review/assets/example.html'):
            self.assertEqual(path.read_bytes(), expected)
        boot = json.loads(expected.decode().split('id="boot" type="application/json">', 1)[1].split('</script>', 1)[0])
        self.assertIsNone(boot['token'])
        self.assertIsNone(boot['paper_demo']['record']['source_utc'])
        self.assertTrue(all(c['teaching'] and not c.get('local_review') and not c.get('reading_events') for c in boot['snapshots']))

    def test_legacy_example_keeps_each_source_response_digest(self):
        html = (ROOT/'docs/legacy-example.html').read_text(encoding='utf-8')
        self.assertEqual(html, (ROOT/'source/brca_review/assets/legacy-example.html').read_text(encoding='utf-8'))
        parts = re.findall(r'响应 SHA-256：<code>([0-9a-f]{64})</code></p><pre[^>]*>([^<]+)</pre>', html)
        self.assertEqual(len(parts), 4)
        for digest, body in parts:
            self.assertEqual(hashlib.sha256(unescape(body).encode()).hexdigest(), digest)

    def test_builder_excludes_unlisted_local_files_and_never_uses_their_timestamps(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            src = base/'source'
            src.mkdir()
            text = b'public example\n'
            (src/'public.txt').write_bytes(text)
            (src/'unlisted-note.txt').write_text('excluded', encoding='utf-8')
            digest = hashlib.sha256(text).hexdigest()
            (src/'MANIFEST.sha256').write_text(digest+'  public.txt\n', encoding='utf-8')
            name, hashed = builder.build(src, base/'a')
            (src/'unlisted-note.txt').write_text('still excluded', encoding='utf-8')
            self.assertEqual(builder.build(src, base/'b'), (name, hashed))
            with zipfile.ZipFile(base/'b'/name) as z:
                self.assertEqual({i.filename for i in z.infolist()}, {builder.STEM+'/'+n for n in ('public.txt','MANIFEST.sha256')})
                self.assertTrue(all(i.date_time == (1980,1,1,0,0,0) and not i.extra and not i.comment for i in z.infolist()))
            self.assertEqual((base/'a'/name).read_bytes(), (base/'b'/name).read_bytes())
            (src/'public.txt').write_bytes(b'changed')
            with self.assertRaises(ValueError):
                builder.build(src, base/'c')
