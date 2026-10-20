"""Self-contained tests for chromatin_potential.py stage A."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from chromatin_potential import (
    CONTACTS,
    intervals_distance,
    point_distance,
)

SCRIPT = Path(__file__).with_name('chromatin_potential.py')
GENE_HEADER = ['gene_name', 'chr', 'start', 'end', 'gene_type', 'from_source', 'N_counts']


def table(path, header, rows):
    with path.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def contact(read, pairtype, start, end, dchr, ds, de, secondary='*', gene='A'):
    return [read, pairtype, 'chr1', start, end, '+', gene, 'lncRNA', 'test', '10M', 0, 60,
            dchr, ds, de, '+', '10M', 0, 30, '*', secondary, '*', '*']


class StageATests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p = Path(self.tmp.name)
        table(self.p/'rna.tsv', GENE_HEADER, [
            ['A', 'chr1', 1, 10000, 'lncRNA', 'test', 10],
            ['B', 'chr1', 1, 10000, 'lncRNA', 'test', 0]])
        table(self.p/'counts.tsv', GENE_HEADER, [
            ['A', 'chr1', 1, 10000, 'lncRNA', 'test', 15],
            ['B', 'chr1', 1, 10000, 'lncRNA', 'test', 0]])
        table(self.p/'uu.tsv', CONTACTS, [
            contact('u1', 'UU', 100, 110, 'chr2', 100, 110),  # trans
            contact('u2', 'UU', 100, 110, 'chr1', 100, 110),  # overlap
            contact('u3', 'UU', 100, 110, 'chr1', 116, 125), # 5 bases gap
            contact('u4', 'UU', 100, 110, 'chr1', 117, 126)]) # 6 bases gap
        table(self.p/'um.tsv', CONTACTS, [
            contact('m1', 'UM', 100, 110, 'chr2', 20, 30, '(chr1,113,+10M,0)'), # pos 114 => gap 3
            contact('m2', 'UM', 100, 110, 'chr2', 20, 30, '(chr1,115,-10M,0)'), # pos 116 => gap 5
            contact('m3', 'UM', 100, 110, 'chr2', 20, 30, '(chr3,0,+10M,0)'),
            contact('m4', 'UM', 100, 110, 'chr1', 100, 110, '(chr3,0,+10M,0)')])

    def run_script(self, mode, *opts):
        return subprocess.run([sys.executable, str(SCRIPT), '--chP_mode', mode,
                               '--rnaseq-counts', str(self.p/'rna.tsv'), '--output-dir', str(self.p/'out'),
                               *map(str, opts)], capture_output=True, text=True)

    def total(self, mode, dist):
        with (self.p/'out'/f'contact_counts_{mode}_dist_{dist}.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(rows), 2)  # includes zero-count gene
        return int(rows[0]['N_counts'])

    def test_boundaries_and_intervals(self):
        self.assertEqual(intervals_distance(100, 110, 116, 125), 5)
        self.assertEqual(point_distance(100, 110, 116), 5)
        self.assertEqual(point_distance(100, 110, 105), 0)
        for dist, expected in ((0, 4), (5, 3), (6, 2), (7, 1)):
            result = self.run_script('unique_dna', '--contacts-uu', self.p/'uu.tsv', '--dist', dist)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(self.total('unique_dna', dist), expected)

    def test_um_filtering(self):
        res = self.run_script('all_dna', '--contacts-uu', self.p/'uu.tsv',
                              '--contacts-um', self.p/'um.tsv', '--dist', 5)
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(self.total('all_dna', 5), 5)  # 3 UU + m2 + m3
        self.assertIn('UM_excluded=2', res.stdout)

    def test_dist_zero_ready_counts_and_fallback(self):
        res = self.run_script('all_dna', '--contact-counts', self.p/'counts.tsv')
        self.assertEqual(res.returncode, 0, res.stderr)

        output = self.p / 'out' / 'chP_all_dna_dist_0.tsv'
        self.assertTrue(output.is_file())

        with output.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))

        self.assertEqual(rows, [])

        res = self.run_script(
            'all_dna',
            '--contacts-uu', self.p/'uu.tsv',
            '--contacts-um', self.p/'um.tsv'
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(self.total('all_dna', 0), 8)

    def test_counts_mode(self):
        res = self.run_script('counts', '--contact-counts', self.p/'counts.tsv')
        self.assertEqual(res.returncode, 0, res.stderr)

        output = self.p / 'out' / 'chP.tsv'
        self.assertTrue(output.is_file())

        with output.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))

        self.assertEqual(rows, [])

    def test_bad_um_and_gene_checks(self):
        with (self.p/'um.tsv').open() as f:
            rows = list(csv.reader(f, delimiter='\t'))
        rows[1][20] = '*'
        table(self.p/'bad.tsv', CONTACTS, rows[1:])
        res = self.run_script('all_dna', '--dist', 1, '--contacts-uu', self.p/'uu.tsv', '--contacts-um', self.p/'bad.tsv')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('bad.tsv:2', res.stderr)
        self.assertIn('requires DNA secondary', res.stderr)
        table(self.p/'extra.tsv', GENE_HEADER, [['UNKNOWN','chr1',1,5,'lncRNA','test',1]])
        res = self.run_script('counts', '--contact-counts', self.p/'extra.tsv')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('genes absent from RNA-seq', res.stderr)

    def test_missing_explicit_file_fails(self):
        res = self.run_script('all_dna', '--contact-counts', self.p/'missing.tsv',
                              '--contacts-uu', self.p/'uu.tsv', '--contacts-um', self.p/'um.tsv')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('file does not exist', res.stderr)


if __name__ == '__main__': unittest.main()
