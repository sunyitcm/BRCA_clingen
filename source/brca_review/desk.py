"""An optional loopback runner for the very same HTML. No third-party dependencies."""
import argparse, datetime as dt, hashlib, json, re, secrets, sys, threading, webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit
from ._engine.brcalite.reading_sources import assignments
from . import ruleprep,erepo,criteria,adoption
from .examples import examples
from .presentation import view_only
from .papers import lookup, stored_xml, EV, PID
from .model import A, BUILDS, LABELS, brief, collect_portal, compare, pack, read
from ._engine.brcalite.__main__ import verify
from ._engine.brcalite.core import parse
from .routes import online
from ._engine.brcalite.reference import Reference
from ._engine.brcalite.workbook import attach_table_notes
from ._engine.brcalite.net import writej

def new_id(): return 's-'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+secrets.token_hex(4)
def render(snaps,index=None,selected=None,token=None):
    snaps=[dict(s,reading_assignments=assignments(s['report']),report=(s['report'] if 'criteria_table' in s['report'] else adoption.effective(s['report'],None))) for s in snaps]
    snaps=[view_only(s) for s in snaps]
    initial=next((s for s in snaps if s['id']==selected),snaps[-1] if snaps else {'id':None,'assembly':'GRCh38'})
    from .guidance import catalog_zh
    boot={'rule_guide':catalog_zh()['genes'],'gene':initial.get('report',{}).get('selected_gene') or initial.get('gene') or 'BRCA1','snapshots':snaps,'index':index if index is not None else list(map(brief,snaps)),'selected':initial['id'],'assembly':initial['assembly'],'labels':LABELS,'token':token,'paper_demo':read(A/'reading_example.json') if any(s.get('teaching') for s in snaps) else {}}
    from html import escape
    rr=initial.get('report') or {};L=rr.get('criteria_table') or {};ep=rr.get('erepo') or {}
    words='；'.join(z['code']+' '+str(z['expert'].get('strength') or '强度未披露') for z in L.get('rows',[]) if (z.get('expert') or {}).get('state')=='Met')
    fallback='<h2>当前选中快照 · '+escape(str(rr.get('gene') or '未选择'))+'</h2><p>'+escape(str((rr.get('facts',{}).get('transcript') or {}).get('hgvsc') or '未核准等位'))+'</p>'+('<p>eRepo 专家补充：'+escape(str(ep['classification']))+'　'+escape(words)+'</p>' if ep.get('state')=='matched' and ep.get('classification') else '')+'<p>本地草案：'+escape(str(rr.get('classification_draft') or '未形成'))+'。启用 JavaScript 可查看全表与本地审核操作。</p>'
    if L.get('rows'):
        subtotal=rr.get('points_accepted')
        fallback+='<p>本基因 '+str(len(L['rows']))+' 条可用主判据；已保存小计：'+escape(str(subtotal if subtotal is not None else '无数值小计'))+'。</p><ul>'+''.join('<li>'+escape(z['code'])+'：'+escape(str((z.get('brief') or {}).get('result') or '暂无取数'))+'</li>' for z in L['rows'] if z['priority'])+'</ul>'
    j=json.dumps(boot,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c').replace('\u2028','\\u2028').replace('\u2029','\\u2029')
    return (A/'shell.html').read_text(encoding='utf-8').replace('@@APP@@',(A/'app.js').read_text(encoding='utf-8').replace('/* REVIEW_LIB */',(A/'criteria.js').read_text(encoding='utf-8'))).replace('@@FALLBACK@@',fallback).replace('@@BOOT@@',j).encode('utf-8')
def check_original(s):
    import base64
    for k,v in s['raw'].items():
        if v['body']:
            b=base64.b64decode(s['blobs'][v['sha256']],validate=True)
            if hashlib.sha256(b).hexdigest()!=v['sha256'] or len(b)!=v['size']:raise ValueError('原文字节校验失败：'+k)
    if hashlib.sha256(base64.b64decode(s['raw_report'],validate=True)).hexdigest()!=s['report_sha256']:raise ValueError('本地报告校验失败')
def build_index(root):
    result=[]
    for f in root.glob('s-*/index.json'):
        try:
            entry=read(f)
            if entry['id']==f.parent.name and (f.parent/'snapshot.json').is_file(): result.append(entry)
        except (KeyError,OSError,ValueError): continue
    return result
class Store:
    def __init__(self,root,workbooks=None,references=None,population_file=None):
        self.dir=root;root.mkdir(parents=True,exist_ok=True);self.workbooks=workbooks;self.population_file=population_file;self.references=references or {};self.demo=examples();ruleprep.configure(root,workbooks)
        self.token=secrets.token_urlsafe(32);self.lock=threading.Lock();self.paper_lock=threading.Lock()
    def with_events(self,s):
        v=dict(s,reading_assignments=assignments(s["report"]));events={}
        for p in (self.dir/'reading_index'/s['id']).glob('*.json'):
            ptr=read(p)
            if re.fullmatch(EV,ptr['event']):
                e=read(self.dir/'reading_events'/ptr['event']/'event.json')
                if e['snapshot']==s['id'] and e['pmid']==p.stem: events[p.stem]=e
        v['reading_events']=events
        return v
    def paper(self,src,pmid,refresh):
        s=self.find(src)
        if not s:raise ValueError('请先选快照')
        e=lookup(self.dir,s,pmid,refresh);folder=self.dir/'reading_events'/e['id']
        s.setdefault('reading_events',{})[pmid]=e
        (folder/'report.html').write_bytes(render([s]))
        idx=self.dir/'reading_index'/s['id'];idx.mkdir(parents=True,exist_ok=True)
        ptr=idx/(e['id']+'.tmp');writej(ptr,{'event':e['id']});ptr.replace(idx/(pmid+'.json'))
        return e
    def paper_bytes(self,i):
        if not re.fullmatch(EV,i):raise ValueError('文献记录不存在')
        e=read(self.dir/'reading_events'/i/'event.json')
        if not e.get('xml_event') or not e.get('article'):raise ValueError('本条没有已核实原文')
        data=stored_xml(self.dir,e['xml_event'],e['pmid']).read_bytes()
        if hashlib.sha256(data).hexdigest()!=e['article']['xml_sha256']:raise ValueError('原文校验失败')
        return data
    def find(self,k,reviewed=True):
        if k is None:return None
        if not isinstance(k,str):raise ValueError('快照 ID 不符')
        for v in self.demo:
            if v['id']==k:return self.project(v,reviewed)
        if not re.fullmatch(r's-[0-9]{8}T[0-9]{6}Z-[a-f0-9]{8}',k):raise ValueError('快照 ID 不符')
        f=self.dir/k/'snapshot.json'
        if not f.is_file():raise ValueError('快照不存在')
        return self.project(read(f),reviewed)
    def checked_base(self,r):
        # Old captures may have F:N but not D:E. Re-read only a known local locked workbook; never fetch during viewing or editing.
        from copy import deepcopy
        if not any((x.get('source') or {}).get('table')=='T9' and not x['source'].get('checked_assignment') for x in r.get('candidates',[])):return r
        path=self.dir/'guidelines/GN092-GN097-v1.2'
        name=next((a['name'] for a in criteria.getj(criteria.R/'rules/attachments.lock.json') if a['gene']==r.get('gene') and a['name'].startswith('Specifications_Table9')),None)
        for d in (path,self.workbooks):
            if d and name and (Path(d)/name).is_file():
                view=deepcopy(r)
                try:attach_table_notes(view,Path(d));return view
                except (ValueError,KeyError,OSError):continue
        return r
    def project(self,s,reviewed):
        base=self.with_events(s)
        if not reviewed:return base
        rev=adoption.latest(self,s['id'])
        if rev and (s['report_sha256']!=rev['base_sha256'] or s['caid']!=rev['caid']):raise ValueError('审核与快照原报告不匹配')
        if rev and rev.get('saved_report_sha256'):
            raw=(adoption.root_for(self,s['id'])/rev['revision']/'report.json').read_bytes()
            if hashlib.sha256(raw).hexdigest()!=rev['saved_report_sha256']:raise ValueError('已保存的审核报告与审核记录哈希不符')
            base['report']=json.loads(raw)
        else:base['report']=adoption.effective(self.checked_base(base['report']),rev)
        if base['report'].get('criteria_table',{}).get('schema')!='brca.criteria.v2':
            from copy import deepcopy
            b=deepcopy(base['report']);b['criteria_table']=criteria.ledger(b,(rev or {}).get('decisions',[]));b['criteria_table']['reviewable']=False;b['legacy_read_only']=True;base['report']=b
        base['local_review']=rev;return view_only(base)
    def review(self,p):
        s=self.find(p['snapshot'],False)
        if not s:raise ValueError('须选定快照')
        check_original(s)
        import base64
        original=json.loads(base64.b64decode(s['raw_report'],validate=True))
        if original.get('caid')!=s['caid'] or original['input']['assembly']!=s['assembly'] or original.get('gene')!=s['gene']:raise ValueError('快照索引与原始报告不一致')
        s=dict(s,report=self.checked_base(original))
        return adoption.save(self,s,p,render)
    def snapshot(self,original,assembly,prior,gene):
        if assembly not in BUILDS or not isinstance(original,str) or not 1<=len(original.strip())<=16384:raise ValueError('请确认组装并输入一个核酸等位')
        v=criteria.scope(gene,assembly,original);old=self.find(prior);sid=new_id();folder=self.dir/sid
        try:wbdir=ruleprep.ensure()
        except (OSError,ValueError):wbdir=None
        reference=Reference(self.references[assembly],assembly) if assembly in self.references else None
        r=online(folder,v,True,reference,workbooks=wbdir is not None)
        r['selected_gene']=gene
        from . import __version__
        r['review_desk_version']=__version__;r['review_policy_sha256']=adoption.policy_hash()
        if r['status']=='review_required' and not criteria.verified(r,gene):
            r['status']='blocked_identity_or_annotation';r['limitations'].append('结果与选择的基因或规范转录本不一致')
        if wbdir is not None and r['status']=='review_required':
            try: attach_table_notes(r,wbdir)
            except (ValueError, OSError, KeyError) as e:
                r['limitations'].append('未能从指定工作簿校验原行：'+str(e))
            writej(folder/'report.json',r)
        if self.population_file and r['status']=='review_required':
            from .selection import population_from_file
            try:population_from_file(self.population_file,r)
            except (ValueError,OSError,KeyError,TypeError) as e:r['source_status']['population_file']='unavailable';r['limitations'].append('本地人群资料未核准，其他证据保留。')
        r['erepo']=erepo.collect(folder,r)
        r.pop('not_applicable_by_spec',None);r['criteria_table']=criteria.ledger(r)
        writej(folder/'report.json',r)
        collect_portal(folder,r)
        s=pack(folder,sid);s['reading_assignments']=assignments(s['report']);s['diff']=compare(old,s);check_original(s)
        writej(folder/'snapshot.json',s)
        (folder/'report.html').write_bytes(render([s]))
        # The index is a commit marker only after the raw data, report and UI exist.
        writej(folder/'index.json',brief(s));return s
def handler(st):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*_):pass
        def send_data(self,code,b,mime,filename=None):
            self.send_response(code);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(b)))
            for k,v in [('Cache-Control','no-store'),('Referrer-Policy','no-referrer'),('X-Content-Type-Options','nosniff'),('X-Frame-Options','DENY')]:self.send_header(k,v)
            if filename:self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
            self.end_headers()
            try:self.wfile.write(b)
            except (BrokenPipeError,ConnectionResetError):pass
        def ans(self,o,status=200):self.send_data(status,json.dumps(o,ensure_ascii=False,separators=(',',':')).encode(),'application/json; charset=utf-8')
        def allowed(self):
            port=self.server.server_port
            return self.headers.get('Host') in (f'127.0.0.1:{port}',f'localhost:{port}')
        def do_GET(self):
            if not self.allowed():return self.ans({'error':'仅限本机'},403)
            path=urlsplit(self.path).path
            if path=='/' and not urlsplit(self.path).query:
                return self.send_data(200,(A/'landing.html').read_bytes(),'text/html; charset=utf-8')
            if path=='/legacy-example':return self.send_data(200,(A/'legacy-example.html').read_bytes(),'text/html; charset=utf-8')
            if path=='/example':return self.send_data(200,(A/'example.html').read_bytes(),'text/html; charset=utf-8')
            if path in ('/','/snapshot','/query'):
                local=sorted(build_index(st.dir),key=lambda x:x['utc']);idx=[brief(v) for v in st.demo]+local
                if path=='/snapshot' and not idx:return self.send_data(200,(A/'example.html').read_bytes(),'text/html; charset=utf-8')
                snaps=[st.find(v['id']) for v in st.demo]
                if local:snaps.append(st.find(local[-1]['id']))
                return self.send_data(200,render(snaps,idx,token=st.token if path!='/snapshot' else None),'text/html; charset=utf-8')
            if path.startswith('/api/paper-xml/'):
                ident=path.partition('/api/paper-xml/')[2]
                try: return self.send_data(200,st.paper_bytes(ident),'application/xml',ident+'.xml')
                except (OSError,KeyError,ValueError):return self.ans({'error':'没有可下载的已校验 XML'},404)
            if path.startswith('/api/snapshot/'):
                try:return self.ans(st.find(path.partition('/api/snapshot/')[2]))
                except ValueError:return self.ans({'error':'找不到快照'},404)
            return self.ans({'error':'没有此页面'},404)
        def do_POST(self):
            if not self.allowed() or not secrets.compare_digest(self.headers.get('X-Review-Token',''),st.token):return self.ans({'error':'请从本机审核页发起'},403)
            origin=self.headers.get('Origin')
            if origin and origin!='http://'+self.headers['Host']:return self.ans({'error':'请求来源不符'},403)
            if self.path not in ('/api/query','/api/paper','/api/review'):return self.ans({'error':'没有此操作'},404)
            if self.headers.get_content_type()!='application/json':return self.ans({'error':'只接受 JSON'},415)
            try:
                n=int(self.headers.get('Content-Length','-1'))
                if not 0<n<=(128*1024 if self.path=='/api/review' else 20000 if self.path=='/api/query' else 2048):raise ValueError('请求长度不符')
                p=json.loads(self.rfile.read(n))
                expected=({'snapshot','pmid','refresh'} if self.path=='/api/paper' else {'gene','assembly','variant','previous'} if self.path=='/api/query' else {'snapshot','base_hash','previous_revision','gene','assembly','caid','schema','selected','note','manual'})
                if self.path=='/api/review' and isinstance(p,dict):p.setdefault('note','');p.setdefault('manual',{})
                if not isinstance(p,dict) or set(p)!=expected or self.path=='/api/review' and p.get('schema')!='brca.selection.v2':raise ValueError('请求字段与操作不符')
            except (ValueError,TypeError,UnicodeError) as e:return self.ans({'error':str(e)},400)
            if self.path=='/api/review':
                try:return self.ans(st.review(p))
                except adoption.ReviewChanged as e:return self.ans({'error':str(e)},409)
                except (ValueError,KeyError,TypeError) as e:return self.ans({'error':str(e)},400)
                except OSError:return self.ans({'error':'审核未能写入；请检查本机存储。'},500)
            if self.path=='/api/paper':
                if not st.paper_lock.acquire(False):return self.ans({'error':'另一篇文献仍在获取，主报告可继续使用。'},409)
                try:return self.ans(st.paper(p['snapshot'],p['pmid'],p['refresh']))
                except (ValueError,TypeError) as e:return self.ans({'error':str(e)},400)
                except Exception as e:return self.ans({'error':'未能完成文献留档；可使用页面的人工阅读链接。'+type(e).__name__},500)
                finally:st.paper_lock.release()
            if not st.lock.acquire(False):return self.ans({'error':'前一次查询仍在进行，请稍后重试。'},409)
            try:return self.ans(st.snapshot(p['variant'],p['assembly'],p['previous'],p['gene']))
            except (ValueError,TypeError) as e:return self.ans({'error':str(e)},400)
            except Exception as e:return self.ans({'error':'未完成留档：'+type(e).__name__+'：'+str(e)},500)
            finally:st.lock.release()
    return Handler
def main():
    p=argparse.ArgumentParser(description='BRCA 基因：一个页面，两种组装。Python 标准库；研究试用。')
    p.add_argument('--port',default=8765,type=int);p.add_argument('--data-dir',type=Path,default=Path.home()/'BRCA-review-runs')
    p.add_argument('--workbooks',type=Path,help='含三份冻结 Excel 原表的目录，可免首次下载')
    for b in ('37','38'):p.add_argument('--ref'+b,type=Path,help='对应组装的 BRCA-lite FASTA JSON 清单')
    p.add_argument('--no-browser',action='store_true');p.add_argument('--export-test',type=Path,help='只生成可双击的离线测试页，不启动服务')
    p.add_argument('--self-test',action='store_true',help='离线核对双组装的身份、版本、原始字节和冻结规则')
    p.add_argument('--prepare-rules',action='store_true',help='仅准备或核验冻结工作簿，不发送任何位点')
    p.add_argument('--population-file',type=Path,help='可选，仅本地：冻结的双队列计数与覆盖 JSON')
    a=p.parse_args();verify()
    if a.prepare_rules:
        ruleprep.configure(a.data_dir,a.workbooks);ruleprep.ensure();print('PASS · 三份工作簿与派生行索引的 SHA-256 均与冻结基线一致。');return
    if a.self_test or a.export_test:
        example=(A/'example.html').read_bytes()
        assert 'CA001182' in example.decode() and b'GRCh37' in example and b'GRCh38' in example
        assert all(route in (A/'landing.html').read_text(encoding='utf-8') for route in ['/query','/snapshot'])
        if a.export_test:
            with a.export_test.open('xb') as f:f.write(example)
        if a.self_test:
            demo=examples()
            assert {(x['gene'],x['assembly']) for x in demo}=={(g,b) for g in criteria.GENES for b in BUILDS}
            for v in demo:
                l=v['report']['criteria_table'];assert l['identity_verified'] and len({z['code'] for z in l['rows']})==19 and sum(z['rule_enabled'] for z in l['rows'])==19
            print('PASS · BRCA1/2 × GRCh37/38、两个入口、19 条本基因可用判据及冻结规则；未请求外站。')
        return
    st=Store(a.data_dir,a.workbooks,{k:v for k,v in [('GRCh37',a.ref37),('GRCh38',a.ref38)] if v},a.population_file)
    try:srv=ThreadingHTTPServer(('127.0.0.1',a.port),handler(st))
    except OSError as e:p.exit(2,'无法在本机监听；可用 --port 0 自动选空闲端口。 '+str(e)+'\n')
    url=f'http://127.0.0.1:{srv.server_port}/';print(url,flush=True)
    if not a.no_browser:webbrowser.open(url)
    try:srv.serve_forever()
    except KeyboardInterrupt:pass
    finally:srv.server_close()
