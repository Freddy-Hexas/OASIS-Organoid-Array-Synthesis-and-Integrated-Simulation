"""Exercise bounded parallel submission and atomic batch validation over HTTP."""
from concurrent.futures import ThreadPoolExecutor
from http.server import ThreadingHTTPServer
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event, Thread
from time import monotonic, sleep
from urllib.error import HTTPError
from urllib.request import Request, urlopen
import json

import gdstk
import server


def main():
    original=(server.INPUT_DIR,server.RUNS_DIR,server.POOL)
    release=Event()
    with TemporaryDirectory(prefix='gds_batch_test_') as root:
        base=Path(root)
        inputs=base/'data';inputs.mkdir()
        runs=base/'runs';runs.mkdir()
        for index in range(3):
            lib=gdstk.Library(name=f'BATCH_{index}')
            cell=lib.new_cell(f'CASE_{index}')
            cell.add(gdstk.rectangle((0,0),(100,100),layer=10))
            ascii_path=base/f'case_{index}.gds'
            lib.write_gds(str(ascii_path))
            (inputs/f'case_{index}.gds').write_bytes(ascii_path.read_bytes())
        assert len(list(inputs.glob('*.gds')))==3,list(inputs.iterdir())
        server.INPUT_DIR=inputs
        server.RUNS_DIR=runs
        server.JOBS.clear()
        server.RECENT_JOBS.clear()
        server.POOL=ThreadPoolExecutor(max_workers=2)
        real_run=server._run_job

        def held_run(job_id,*_args):
            server._update_job(job_id,status='running',stage='测试占用工作线程',progress=.2)
            if not release.wait(12):
                raise TimeoutError('Test worker was not released')
            server._update_job(job_id,status='complete',stage='完成',progress=1.,
                               finished_at=server.utc_now(),message='测试完成',
                               result={'capacity_interval':{'lower_bound':1}})

        server._run_job=held_run
        http=ThreadingHTTPServer(('127.0.0.1',0),server.Handler)
        Thread(target=http.serve_forever,daemon=True).start()
        url=f'http://127.0.0.1:{http.server_port}'

        def call(endpoint,payload=None):
            data=None if payload is None else json.dumps(payload).encode('utf-8')
            request=Request(url+endpoint,data=data,
                            headers={'Content-Type':'application/json'} if data else {})
            try:
                with urlopen(request,timeout=30) as response:
                    raw=response.read()
                    return json.loads(raw) if response.headers.get_content_type()=='application/json' else raw
            except HTTPError as exc:
                return exc.code,json.loads(exc.read())

        try:
            files=call('/api/files')['files']
            assert len(files)==3,files
            ids=[item['id'] for item in files]
            assert call('/admin')[:15]==b'<!doctype html>'
            assert call('/admin.js')[:5]==b'const'
            invalid=call('/api/jobs/batch',{'file_ids':ids,'method':'four_side_pads',
                                            'rules':{'wire_width_um':0}})
            assert invalid[0]==422 and not server.JOBS,invalid
            invalid=call('/api/jobs/batch',{'file_ids':ids,
                                            'support_layers':{ids[1]:'999/0'}})
            assert invalid[0]==422 and not server.JOBS,invalid
            submitted=call('/api/jobs/batch',{'file_ids':ids,'method':'four_side_pads',
                                              'rules':{'electrode_diameter_um':20,'wire_width_um':2,
                                                       'spacing_um':0,'margin_um':0,
                                                       'minimum_center_spacing_um':0},
                                              'pad_settings':{'pad_width_um':500}})
            assert submitted['accepted']==3,submitted
            for entry in submitted['jobs']:
                stored=server.JOBS[entry['job_id']]['rules']
                assert stored['electrode_diameter_um']==20 and stored['wire_width_um']==2,stored
                assert stored['spacing_um']==stored['margin_um']==stored['minimum_center_spacing_um']==0,stored
            assert {server.JOBS[j['job_id']]['batch_id'] for j in submitted['jobs']}=={submitted['batch_id']}
            deadline=monotonic()+8
            while monotonic()<deadline:
                queue=call('/api/queue')
                if queue['running']==2 and queue['queued']==1:
                    break
                sleep(.05)
            assert queue['running']==2 and queue['queued']==1,queue
            assert set(queue['jobs'])==set(ids)
            assert call('/api/jobs',{'file_id':ids[0]})[0]==409
            assert call('/api/jobs/batch',{'file_ids':[ids[2],ids[0]]})[0]==409
            assert len(server.JOBS)==3
            release.set()
            deadline=monotonic()+8
            while monotonic()<deadline:
                queue=call('/api/queue')
                if not queue['running'] and not queue['queued']:
                    break
                sleep(.05)
            assert queue['running']==queue['queued']==0,queue
            assert all(job['lower_bound']==1 for job in queue['jobs'].values()),queue
            assert all(call('/api/job/'+entry['job_id'])['status']=='complete'
                       for entry in submitted['jobs'])
            print('Parallel workers, queued third task, duplicate guard, batch validation and admin assets passed',flush=True)
        finally:
            release.set()
            http.shutdown();http.server_close()
            server.POOL.shutdown(wait=True)
            server._run_job=real_run
            server.INPUT_DIR,server.RUNS_DIR,server.POOL=original
            server.JOBS.clear()
            server.RECENT_JOBS.clear()


if __name__=='__main__':
    main()
