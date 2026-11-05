# Checks for the site revision items (B1 to B3) against the local study data.
# Run from the repository root: python -m unittest tests.test_site_revision
# Needs secret.cfg, static/study_data with at least one VDJbase dataset, and network access to MADC at startup.

import os
import sys
import unittest
from glob import glob

sys.path.insert(0, os.getcwd())
from sqlalchemy import func
import app as backend
from db.vdjbase_airr_model import Study, Sample
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH

SPECIES, DATASET = 'Human', 'IGH'
CURL = 'curl -w "HTTP code: %{response_code}" '


def script_lines(resp):
    return [line for line in resp.data.decode('utf-8').splitlines() if line.startswith('curl')]


def local_path(url):
    """The test client path for a URL built from BACKEND_LINK"""
    return '/' + url[len(backend.app.config['BACKEND_LINK']):]


@unittest.skipUnless(SPECIES in backend.vdjbase_dbs and DATASET in backend.vdjbase_dbs[SPECIES], 'no local Human IGH data')
class StudyScriptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.client = backend.app.test_client()
        session = backend.vdjbase_dbs[SPECIES][DATASET].get_session()
        # the smallest study, by sample count
        cls.study = session.query(Study.study_name).join(Sample, Sample.study_id == Study.id) \
            .group_by(Study.study_name).order_by(func.count(Sample.id)).first()[0]

    def get(self, path):
        return self.client.get('/api/repseq/' + path)

    def check_script(self, kind):
        resp = self.get(f'download_study_{kind}/{SPECIES}/{self.study}/{DATASET}')
        self.assertEqual(resp.status_code, 200)
        self.assertIn(f'{self.study}_{DATASET}_{kind}_download.sh', resp.headers['Content-Disposition'])
        lines = script_lines(resp)
        self.assertTrue(lines)
        for line in lines:
            self.assertTrue(line.startswith(CURL))
            url, file_name = line[len(CURL):].split(' > ')
            self.assertEqual(os.path.basename(url), file_name)
            with self.client.get(local_path(url)) as file_resp:
                self.assertEqual(file_resp.status_code, 200, url)
        return lines

    def test_genotype_script(self):
        # one line per sample whose genotype file is on disk; samples without the file are skipped
        on_disk = glob(os.path.join(VDJBASE_SAMPLE_PATH, SPECIES, DATASET, self.study, '*', '*_genotype.tsv'))
        self.assertEqual(len(self.check_script('genotypes')), len(on_disk))

    def test_haplotype_script(self):
        from api.vdjbase.vdjbase import study_haplotype_files, existing_sample_files
        session = backend.vdjbase_dbs[SPECIES][DATASET].get_session()
        if not existing_sample_files(SPECIES, DATASET, study_haplotype_files(session, self.study)):
            self.skipTest(f'{self.study} has no haplotype files')
        self.check_script('haplotypes')

    def test_airr_script_format(self):
        taxid = backend.vdjbase_dbs[SPECIES][DATASET].taxid
        recs = backend.madc_index.get(taxid, {}).get(DATASET, {})
        if not recs:
            self.skipTest('MADC index not loaded')
        study_id = next(iter(recs.values()))['study_id']
        resp = self.get(f'download_study_script/{SPECIES}/{study_id}/{DATASET}')
        self.assertEqual(resp.status_code, 200)
        lines = script_lines(resp)
        self.assertEqual(len(lines), sum(1 for r in recs.values() if r['study_id'] == study_id))
        self.assertTrue(all(line.startswith(CURL + 'https://madc.vdjbase.org/airr/v1/rearrangement/') for line in lines))

    def test_unknown_returns_404(self):
        for kind in ('script', 'genotypes', 'haplotypes'):
            self.assertEqual(self.get(f'download_study_{kind}/Nope/{self.study}/{DATASET}').status_code, 404)
            self.assertEqual(self.get(f'download_study_{kind}/{SPECIES}/{self.study}/Nope').status_code, 404)
        for kind in ('genotypes', 'haplotypes'):
            self.assertEqual(self.get(f'download_study_{kind}/{SPECIES}/NoSuchStudy/{DATASET}').status_code, 404)


if __name__ == '__main__':
    unittest.main()
