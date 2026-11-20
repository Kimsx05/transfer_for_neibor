"""Read-only, on-demand individual cards for every source-linked cell and unlocalized region."""
from cards import Cards
from base6 import *
from http.server import BaseHTTPRequestHandler,HTTPServer
from urllib.parse import urlparse,parse_qs,quote
import io,html
if __name__=='__main__':
 section=sys.argv[1];c=Cards(section);keys=list(c.q.index)
 class Handler(BaseHTTPRequestHandler):
  def do_GET(self):
   url=urlparse(self.path);q=parse_qs(url.query);i=int(q.get('i',['0'])[0]);i=max(0,min(len(keys)-1,i));k=q.get('k',[keys[i]])[0]
   if k not in c.q.index:self.send_error(404);return
   if url.path=='/card':
    buf=io.BytesIO();c.card(k).save(buf,format='JPEG',quality=92);body=buf.getvalue();content='image/jpeg'
   elif url.path=='/':
    i=c.q.index.get_loc(k);body=f'<meta charset=utf-8><title>{section} | {html.escape(k)}</title><p>{section}, {i+1} / {len(keys)}. Nucleus contours, inferred assignment, observed grid counts.</p><a href="/?i={max(0,i-1)}">Previous</a> · <a href="/?i={min(len(keys)-1,i+1)}">Next</a><form>Review ID <input name=k value="{html.escape(k)}"><button>Open</button></form><img src="/card?k={quote(k)}">'.encode();content='text/html; charset=utf-8'
   else:self.send_error(404);return
   self.send_response(200);self.send_header('Content-Type',content);self.send_header('Content-Length',str(len(body)));self.end_headers();self.wfile.write(body)
 HTTPServer(('127.0.0.1',int(sys.argv[2]) if len(sys.argv)>2 else 8765),Handler).serve_forever()
