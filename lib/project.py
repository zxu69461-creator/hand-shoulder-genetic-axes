"""Results database shared by the analysis scripts.

Each script writes its result tables to one DuckDB file (config.DB); downstream scripts and the Source Data read them
from there. Large tables are written as Parquet next to the database and exposed as views.

    P = Project(__file__)
    with P.run(__file__, params={...}) as run:
        run.input('GCST002783', path=..., note='...')     # records an input (printed)
        run.table('mr_results', df, note='...')            # CREATE OR REPLACE TABLE
        run.parquet('factor_gwas', df, note='...')         # Parquet file + view
    with P.connect() as con:
        con.execute('select ...').df()
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402


class Run:
    def __init__(self, project, script, params):
        self.p, self.script, self.params = project, Path(script).name, params or {}

    def __enter__(self):
        print(f'== {self.script}', {k: v for k, v in self.params.items() if not isinstance(v, (dict, list))} or '', flush=True)
        return self

    def __exit__(self, *exc):
        return False

    def input(self, name, path=None, url=None, note=None):
        print(f'   input {name}' + (f' ({note})' if note else ''), flush=True)

    @staticmethod
    def _drop(con, name):
        """remove an earlier table or view of that name"""
        if con.execute('SELECT count(*) FROM duckdb_views() WHERE view_name = ?', [name]).fetchone()[0]:
            con.execute(f'DROP VIEW {name}')
        con.execute(f'DROP TABLE IF EXISTS {name}')

    def table(self, name, df, note=None):
        with self.p.connect() as con:
            con.register('_df', df)
            self._drop(con, name)
            con.execute(f'CREATE TABLE {name} AS SELECT * FROM _df')
            con.unregister('_df')
        print(f'   table {name}: {len(df):,} rows', flush=True)
        return len(df)

    def parquet(self, name, df, note=None):
        d = self.p.db.parent / 'parquet'
        d.mkdir(parents=True, exist_ok=True)
        f = d / f'{name}.parquet'
        with self.p.connect() as con:
            con.register('_df', df)
            con.execute(f"COPY (SELECT * FROM _df) TO '{f}' (FORMAT PARQUET, COMPRESSION ZSTD)")
            con.unregister('_df')
            self._drop(con, name)
            con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{f}')")
        print(f'   parquet {name}: {len(df):,} rows', flush=True)
        return len(df)


class Project:
    def __init__(self, script):
        self.root = Path(script).resolve().parents[1]
        self.db = C.DB

    def connect(self, read_only=False):
        """Connection to the results database (one configuration per process, so `read_only` is accepted but not used)."""
        self.db.parent.mkdir(parents=True, exist_ok=True)
        return duckdb.connect(str(self.db))

    def run(self, script, params=None):
        return Run(self, script, params)


def fetch(url, dest):
    """Download url to dest with curl (resumable); returns dest."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        subprocess.run(['curl', '-fL', '--retry', '5', '-C', '-', '-o', str(dest) + '.part', url], check=True)
        Path(str(dest) + '.part').replace(dest)
    return dest
