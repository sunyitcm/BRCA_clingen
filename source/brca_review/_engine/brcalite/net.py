"""Bounded read-only public HTTP; record success, partial failures and exact bytes."""
from xml.etree import ElementTree as ET
import datetime as dt, hashlib, json, time, urllib.error, urllib.request, socket, http.client
from pathlib import Path
from urllib.parse import urlsplit
from email.utils import parsedate_to_datetime
import math

_host_cooldown = {}

def retry_seconds(value, now=None):
    """Seconds until a service permits retry, including fractional vendor values."""
    if value is None:
        return None
    try:
        seconds = float(value)
        return seconds if math.isfinite(seconds) and seconds >= 0 else None
    except (ValueError, TypeError):
        try:
            when = parsedate_to_datetime(value)
            if when.tzinfo is None:
                return None
            return max(0.0, when.timestamp() - (time.time() if now is None else now))
        except (ValueError, TypeError, OverflowError):
            return None

class SameServiceRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        old, new = urlsplit(req.full_url), urlsplit(newurl)
        if new.scheme != 'https' or new.netloc != old.netloc:
            raise ValueError('unexpected redirect to another service')
        return super().redirect_request(req, fp, code, msg, headers, newurl)

def sha(x): return hashlib.sha256(x).hexdigest()
def utc(): return dt.datetime.now(dt.timezone.utc).isoformat()
def writej(p,o): Path(p).write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def getj(p): return json.loads(Path(p).read_bytes())
def receipt_headers(headers):
    keep={'content-type','content-length','date','last-modified','etag','retry-after','content-encoding','x-ratelimit-limit','x-ratelimit-remaining','x-ratelimit-reset'}
    return [(k,v) for k,v in headers.items() if k.lower() in keep]
def fetch(root,label,url,payload=None,kind='json',attempts=2,timeout=15):
    target=Path(root)/label; target.mkdir(parents=True,exist_ok=False)
    data=None if payload is None else json.dumps(payload,separators=(',',':')).encode()
    from . import __version__
    headers={'Accept':{'json':'application/json','xml':'application/xml','text':'text/plain','html':'text/html'}[kind],'User-Agent':f'brca-lite/{__version__} research (public read-only)'}
    if data is not None: headers['Content-Type']='application/json'
    meta={'label':label,'method':'GET' if data is None else 'POST','url':url,'request_headers':headers,'request_sha256':None if data is None else sha(data),'at_utc':utc(),'attempts':[]}
    if data is not None:(target/'request.json').write_bytes(data)
    host=urlsplit(url).netloc
    started=time.monotonic()
    for n in range(attempts):
        delay=_host_cooldown.get(host,0)-time.time()
        if delay>12:
            meta.update(ok=False,state='host_cooldown',retry_at_utc=dt.datetime.fromtimestamp(_host_cooldown[host],dt.timezone.utc).isoformat(),finished_utc=utc(),elapsed_s=round(time.monotonic()-started,4))
            writej(target/'receipt.json',meta)
            return None
        if delay>0:time.sleep(delay)
        current={'at_utc':utc()}; raw=None; status=None; retry_after=None; network_start=time.monotonic()
        transient=False
        try:
            req=urllib.request.Request(url,data=data,headers=headers)
            with urllib.request.build_opener(SameServiceRedirect()).open(req,timeout=timeout) as r:
                status=r.status; current.update(status=status,final_url=r.geturl(),response_headers=receipt_headers(r.headers))
                raw=r.read(12*1024*1024+1)
                if len(raw)>12*1024*1024:raise ValueError('response exceeds limit')
                if kind=='json':
                    obj=json.loads(raw)
                    if isinstance(obj,dict) and (obj.get('error') or obj.get('errorType') or obj.get('errorCode')):raise ValueError('application error')
                elif kind=='xml':
                    elem=ET.fromstring(raw)
                    if elem.tag.lower() in ('html','error') or elem.find('.//ERROR') is not None:raise ValueError('not a valid XML service response')
                current['ok']=True
        except urllib.error.HTTPError as e:
            status=e.code;current.update(status=status,error=str(e),response_headers=receipt_headers(e.headers))
            try:
                raw=e.read(131072)
            except (OSError,http.client.HTTPException) as read_error:
                current['error_body_read']=str(read_error)
            finally:
                e.close()
            transient=status==429 or 500<=status<=599
            retry_after=retry_seconds(e.headers.get('Retry-After'))
        except (OSError,http.client.HTTPException) as e:
            current.update(error=str(e),error_kind='transport');transient=True
            if isinstance(e,http.client.IncompleteRead):raw=e.partial
        except (ValueError,ET.ParseError) as e:current['error']=str(e)
        current.update(elapsed_s=round(time.monotonic()-network_start,4),finished_utc=utc())
        wait=max(2**n,retry_after or 0)
        if transient and (status==429 or retry_after is not None):
            _host_cooldown[host]=time.time()+wait
            current['retry_at_utc']=dt.datetime.fromtimestamp(_host_cooldown[host],dt.timezone.utc).isoformat()
        if raw is not None:
            name='body' if current.get('ok') else f'failed-{n+1}.body'
            (target/name).write_bytes(raw);current.update(body_file=name,body_sha256=sha(raw),bytes=len(raw))
        meta['attempts'].append(current);meta['ok']=bool(current.get('ok'));meta['finished_utc']=utc()
        meta['elapsed_s']=round(time.monotonic()-started,4)
        writej(target/'receipt.json',meta)
        if meta['ok']:return raw
        if not transient or wait>12:break
        if n+1<attempts:time.sleep(wait)
    return None
