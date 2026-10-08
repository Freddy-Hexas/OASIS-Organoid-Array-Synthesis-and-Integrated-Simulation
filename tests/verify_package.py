"""Check a relocated package, HTTP assets, and a real electrode-to-Pad job.

Run from the package root: python -B tests/verify_package.py [--ui]
All servers, generated inputs, caches and results use temporary directories.
No production case is routed. No solver time limit is set.
"""
from argparse import ArgumentParser
from contextlib import contextmanager
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from time import monotonic, sleep
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import gzip
import json
import shutil
import socket
import subprocess
import sys
import os

import gdstk


PACKAGE = Path(__file__).resolve().parents[1]


def request(base, endpoint, payload=None):
    data = None if payload is None else json.dumps(payload).encode('utf-8')
    headers = {'Content-Type': 'application/json'} if data is not None else {}
    with urlopen(Request(base + endpoint, data=data, headers=headers), timeout=60) as response:
        raw = response.read()
        if response.headers.get('Content-Encoding') == 'gzip':
            raw = gzip.decompress(raw)
        if response.headers.get_content_type() == 'application/json':
            return json.loads(raw)
        return raw


@contextmanager
def serve(package, inputs, runs, workdir, logdir):
    with socket.socket() as probe:
        probe.bind(('127.0.0.1', 0))
        port = probe.getsockname()[1]
    url = f'http://127.0.0.1:{port}'
    stdout = logdir / f'{port}.stdout.log'
    stderr = logdir / f'{port}.stderr.log'
    env = dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1')
    with stdout.open('w', encoding='utf-8') as out, stderr.open('w', encoding='utf-8') as err:
        child = subprocess.Popen([sys.executable, '-B', str(package / 'start_workbench.py'),
            '--port', str(port), '--workers', '2', '--input-dir', str(inputs),
            '--output-dir', str(runs)], cwd=workdir, env=env, stdout=out, stderr=err)
        try:
            # This readiness limit checks server startup only, not routing.
            until = monotonic() + 30
            while True:
                if child.poll() is not None:
                    raise AssertionError(stderr.read_text(encoding='utf-8'))
                try:
                    health = request(url, '/api/health')
                    if health['ok']:
                        break
                except (URLError, TimeoutError):
                    if monotonic() >= until:
                        raise AssertionError('Server did not start: ' + stderr.read_text(encoding='utf-8'))
                    sleep(.1)
            assert Path(health['input_root']) == inputs.resolve()
            yield url
        finally:
            child.terminate()
            try:
                child.wait(timeout=20)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()


def browser_check(url, count, *, result_count=None, screenshot=None):
    from playwright.sync_api import sync_playwright, expect
    errors = []
    with sync_playwright() as runtime:
        browser = runtime.chromium.launch(headless=True)
        page = browser.new_page(viewport={'width': 1440, 'height': 1000}, device_scale_factor=1)
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(url, wait_until='domcontentloaded')
        expect(page.locator('#fileList button')).to_have_count(count, timeout=60000)
        expect(page.locator('#runButton')).to_be_enabled(timeout=60000)
        expect(page.locator('#vectorViewport canvas')).to_be_visible(timeout=60000)
        if result_count is None:
            expect(page.locator('#minimumCenterSpacingMm')).to_have_value('0.07')
            expect(page.locator('#electrodeRegionRadiusMm')).to_have_value('3')
        else:
            expect(page.locator('#electrodeHeroCount')).to_have_text(str(result_count), timeout=60000)
            expect(page.locator('#stageCards button')).to_have_count(5)
            expect(page.locator('#distributionCount')).to_have_text(str(result_count), timeout=60000)
            page.locator('[data-formula="maximum"]').first.click()
            expect(page.locator('#distributionFormulaDialog')).to_be_visible()
            assert page.locator('#distributionFormulaBody .katex').count() > 0
            assert page.locator('#distributionFormulaBody .katex-error').count() == 0
            page.locator('#distributionFormulaClose').click()
        # Exercise the actual canvas pointer path used for zoom and pan.
        canvas = page.locator('#vectorViewport canvas')
        initial = canvas.screenshot()
        box = canvas.bounding_box()
        x, y = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
        page.mouse.move(x, y)
        page.mouse.wheel(0, -300)
        page.wait_for_timeout(300)
        page.mouse.down()
        page.mouse.move(x + 70, y + 40, steps=8)
        page.mouse.up()
        page.wait_for_timeout(300)
        assert canvas.screenshot() != initial, 'Canvas did not respond to zoom/pan'
        if screenshot:
            page.screenshot(path=str(screenshot), full_page=True)
        page.goto(url + '/admin', wait_until='domcontentloaded')
        expect(page.locator('#adminFileList input[type="checkbox"]')).to_have_count(count, timeout=60000)
        page.locator('#selectAll').click()
        expect(page.locator('#submitBatch')).to_be_enabled()
        page.locator('#clearAll').click()
        expect(page.locator('#submitBatch')).to_be_disabled()
        browser.close()
    assert not errors, errors
    return {'page_errors': errors, 'canvas_zoom_pan': True, 'admin_selection': True,
            'formula_rendered': result_count is not None}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--ui', action='store_true', help='Also check Chromium UI (requires Playwright)')
    parser.add_argument('--output', type=Path, help='Optional verification report, outside the temporary workspace')
    args = parser.parse_args()
    report = {'passed': False, 'production_cases_routed': 0, 'solver_time_limit_added': False}
    with TemporaryDirectory(prefix='gds_portable_package_') as temp:
        base = Path(temp)
        package = base / 'portable-package'
        shutil.copytree(PACKAGE, package,
                        ignore=shutil.ignore_patterns('outputs', '__pycache__', '.venv', '.git', '*.pyc'))
        different_cwd = base / 'unrelated working directory'
        different_cwd.mkdir()
        data = package / 'data'
        runs = base / 'sample runs'
        names = sorted(p.name for p in data.iterdir() if p.is_file() and p.suffix.lower() == '.gds')
        assert names, 'No packaged input structures'
        print('Checking relocated startup, all input GDS files and browser assets', flush=True)
        with serve(package, data, runs, different_cwd, base) as url:
            records = request(url, '/api/files')
            assert sorted(p['name'] for p in records['files']) == names
            assert records['tasks'] == {}
            assert request(url, '/api/queue')['max_workers'] == 2
            checked = []
            for item in records['files']:
                assert item['relative_path'] == 'data/' + item['name']
                assert item['sha256'] == sha256((data / item['name']).read_bytes()).hexdigest()
                detail = request(url, '/api/file/' + item['id'])
                assert detail['suggested_support_layer'] in detail['layers']
                assert detail['reference'] is None
                checked.append({'name': item['name'], 'sha256': item['sha256'],
                                'support_layer': detail['suggested_support_layer']})
            asset_root = package / 'gds_frontend' / 'web_app'
            assets = {'/': 'index.html', '/admin': 'admin.html',
                      **{('/' + name): name for name in ('app.js', 'app.css', 'admin.js', 'admin.css',
                         'distribution.js', 'distribution.css', 'vector-viewer.js', 'vector-worker.js',
                         'vendor/katex/katex.min.js', 'vendor/katex/katex.min.css',
                         'vendor/katex/fonts/KaTeX_Main-Regular.woff2')}}
            for endpoint, local in assets.items():
                assert request(url, endpoint) == (asset_root / local).read_bytes(), endpoint
            example = records['files'][0]
            scene = request(url, '/api/scene/input/' + example['id'])
            assert scene['pathCount'] > 0 and len(scene['paths']) == scene['pathCount']
            assert b'<svg' in request(url, '/api/vector/input/' + example['id'])
            report['inputs'] = checked
            report['http_assets_checked'] = len(assets)
            report['different_working_directory'] = True
            report['relocated_package'] = True
            if args.ui:
                report['input_ui'] = browser_check(url, len(names))
        # This unfamiliar input has no generator file or preset connectivity.
        custom = base / 'custom inputs'
        custom.mkdir()
        source = custom / 'unfamiliar_strip.GDS'
        library = gdstk.Library(unit=1e-6, precision=1e-9)
        cell = library.new_cell('UNLABELED_SUPPORT')
        cell.add(gdstk.rectangle((534, -452), (1214, -412), layer=42),
                 gdstk.rectangle((1254, -452), (1934, -412), layer=42))
        library.write_gds(str(source))
        custom_runs = base / 'custom results'
        print('Running one synthetic input through placement, Pad routing and exported-GDS audit', flush=True)
        with serve(package, custom, custom_runs, different_cwd, base) as url:
            item = request(url, '/api/files')['files'][0]
            assert item['relative_path'] == 'data/unfamiliar_strip.GDS'
            submitted = request(url, '/api/jobs', {'file_id': item['id'],
                'support_layer': '42/0', 'method': 'four_side_pads', 'rules': {
                'electrode_region_radius_um': 40, 'minimum_center_spacing_um': 120, 'margin_um': 1}})
            job_id = submitted['job_id']
            previous_stage = None
            while True:
                job = request(url, '/api/job/' + job_id)
                if job.get('stage') != previous_stage:
                    print('Stage:', job.get('stage'), flush=True)
                    previous_stage = job.get('stage')
                if job['status'] not in ('queued', 'running'):
                    break
                sleep(.5)
            assert job['status'] == 'complete', job.get('message')
            routing = job['result']['routing']
            assert routing['retained_routes'] == 1
            assert routing['gds_roundtrip_audit']['passed']
            assert routing['integer_polygon_audit']['passed']
            assert job['result']['capacity_interval']['lower_bound'] == 1
            assert routing['settings']['milp_time_limit_s'] is None
            for artifact in job['artifacts']:
                raw = request(url, '/api/artifact/' + job_id + '/' + artifact)
                if isinstance(raw, (dict, list)):
                    assert raw == json.loads((custom_runs / job_id / artifact).read_text(encoding='utf-8'))
                else:
                    assert raw == (custom_runs / job_id / artifact).read_bytes()
            for kind in ('job', 'feasible'):
                assert request(url, '/api/scene/' + kind + '/' + job_id)['pathCount'] > 0
            has_curves = any(route.get('curve', {}).get('curve_segments') for route in routing['routes'])
            if has_curves:
                assert request(url, '/api/scene/curve/' + job_id)['pathCount'] > 0
            else:
                try:
                    request(url, '/api/scene/curve/' + job_id)
                except HTTPError as error:
                    assert error.code == 404
                else:
                    raise AssertionError('A straight route exposed a curved-centerline view')
            preview = request(url, '/api/job/' + job_id + '?detail=preview')
            assert preview['detail_level'].startswith('preview')
            metrics = request(url, '/api/distribution/' + job_id +
                '?radius_um=100&center_x_um=1234&center_y_um=-432&resolution=192')
            assert metrics['electrode_count'] == 1 and metrics['heatmap']
            report['synthetic_job'] = {'connected_electrodes': 1, 'roundtrip_audit': True,
                'integer_polygon_audit': True, 'distribution_api': True,
                'vector_stages': True, 'curve_view_matches_geometry': True,
                'artifacts_downloaded': len(job['artifacts'])}
            if args.ui:
                screenshot = None
                if args.output:
                    args.output.parent.mkdir(parents=True, exist_ok=True)
                    screenshot = args.output.with_suffix('.png')
                report['result_ui'] = browser_check(url, 1, result_count=1, screenshot=screenshot)
        with serve(package, custom, custom_runs, different_cwd, base) as url:
            restored = request(url, '/api/files')['tasks'][item['id']]
            assert restored['id'] == job_id and restored['status'] == 'complete' and restored['lower_bound'] == 1
            assert request(url, '/api/job/' + job_id)['status'] == 'complete'
            report['restart_restores_results'] = True
        assert not (package / 'gds_frontend' / 'web_app' / 'runs').exists()
        assert not (package / 'gds_frontend' / 'results').exists()
        report['legacy_output_paths_absent'] = True
    report['passed'] = True
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
