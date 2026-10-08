"""Start the GDS workbench from any working directory.

The launcher configures storage and HTTP settings only. It does not change
the process rules or impose a solver time limit.
"""
from argparse import ArgumentParser
from pathlib import Path
import os
import runpy


PACKAGE_DIR = Path(__file__).resolve().parent


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--host', default='127.0.0.1', help='Bind address (default: local computer only)')
    parser.add_argument('--port', type=int, default=8769, help='HTTP port (default: 8769)')
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=2,
                        help='Parallel jobs, 1 to 4 (default: 2)')
    parser.add_argument('--input-dir', type=Path, default=PACKAGE_DIR / 'data',
                        help='Folder of input GDS files (default: package/data)')
    parser.add_argument('--output-dir', type=Path, default=PACKAGE_DIR / 'outputs' / 'runs',
                        help='Run and cache folder (default: package/outputs/runs)')
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error('--port must be between 1 and 65535')
    if not args.host.strip():
        parser.error('--host must not be empty')
    inputs = args.input_dir.expanduser().resolve()
    outputs = args.output_dir.expanduser().resolve()
    if not inputs.is_dir():
        parser.error(f'Input folder does not exist: {inputs}')
    if inputs == outputs:
        parser.error('Input and output folders must be different')
    outputs.mkdir(parents=True, exist_ok=True)
    settings = {'GDS_WORKBENCH_HOST': args.host.strip(),
                'GDS_WORKBENCH_PORT': str(args.port),
                'GDS_WORKBENCH_WORKERS': str(args.workers),
                'GDS_WORKBENCH_DATA_DIR': str(inputs),
                'GDS_WORKBENCH_RUNS_DIR': str(outputs)}
    os.environ.update(settings)
    runpy.run_path(str(PACKAGE_DIR / 'gds_frontend' / 'web_app' / 'server.py'), run_name='__main__')


if __name__ == '__main__':
    main()
