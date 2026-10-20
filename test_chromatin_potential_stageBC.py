"""Regression tests for chromatin_potential.py stages B and C.

Run next to chromatin_potential.py:
    python -m unittest -v test_chromatin_potential_stageBC.py
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest


SCRIPT = Path(__file__).with_name('chromatin_potential.py')
HEADER = ('gene_name', 'chr', 'start', 'end', 'gene_type', 'from_source', 'N_counts')
OUTPUT_COLUMNS = ('gene_name', 'gene_type', 'N_counts_RNAseq_plus_1',
                  'gene_length', 'N_contacts_plus_1', 'chP', 'pval', 'fdr_bh')
CONTACT_COLUMNS = ('read_id', 'pairtype', 'rna_chr', 'rna_start', 'rna_end',
                   'rna_strand', 'gene_name', 'gene_type', 'from_source',
                   'rna_cigar', 'rna_NM', 'rna_mapq', 'dna_chr', 'dna_start',
                   'dna_end', 'dna_strand', 'dna_cigar', 'dna_NM', 'dna_mapq',
                   'rna_secondary_alignments', 'dna_secondary_alignments',
                   'rna_other_tags', 'dna_other_tags')


def write_tsv(path, columns, rows):
    """Write a small input fixture."""
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(rows)


def read_tsv(path):
    """Read a TSV result, including its header."""
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        return tuple(reader.fieldnames), list(reader)


def contact(read_id, gene, biotype, dna_start, pairtype='UU', secondary='*'):
    """Create a valid Voting contact record."""
    return [read_id, pairtype, 'chr1', 100, 109, '+', gene, biotype,
            'test', '10M', 0, 60, 'chr1', dna_start, dna_start + 9,
            '+', '10M', 0, 60, '*', secondary, '*', '*']


class StageBCTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.rna = self.root / 'rna.tsv'
        self.contacts = self.root / 'counts.tsv'
        self.uu = self.root / 'uu.tsv'
        self.out = self.root / 'out'

        # Gene D is a zero-contact RNA-seq anchor. Gene C fails length > 200.
        write_tsv(self.rna, HEADER, [
            ['A', 'chr1', 1, 1000, 'protein_coding', 'test', 10],
            ['B', 'chr1', 1, 1000, 'lncRNA', 'test', 150],
            ['C', 'chr1', 1, 200, 'protein_coding', 'test', 3],
            ['D', 'chr1', 1, 1000, 'lncRNA', 'test', 0],
        ])
        write_tsv(self.contacts, HEADER, [
            ['A', 'chr1', 1, 1000, 'protein_coding', 'test', 150],
            ['B', 'chr1', 1, 1000, 'lncRNA', 'test', 20],
            ['C', 'chr1', 1, 200, 'protein_coding', 'test', 300],
        ])

    def run_script(self, mode='counts', *options, plot=False, out=None):
        """Execute the actual CLI, not a copy of its implementation."""
        args = [sys.executable, str(SCRIPT), '--chP_mode', mode,
                '--rnaseq-counts', str(self.rna), '--output-dir', str(out or self.out),
                '--N_contacts_min', '0', '--gene_len_min', '200',
                '--fdr_threshold', '0.05']
        if mode == 'counts':
            args.extend(['--contact-counts', str(self.contacts)])
        if plot:
            args.extend(['--protein-coding-biotype', 'protein_coding'])
        return subprocess.run(args + list(map(str, options)),
                              capture_output=True, text=True, timeout=90)

    def test_z_statistics_and_bh_fdr(self):
        """Verify exact Z, pseudocounts, gene lengths, and BH adjustment."""
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        columns, output = read_tsv(self.out / 'chP.tsv')
        self.assertEqual(columns, OUTPUT_COLUMNS)
        self.assertEqual({row['gene_name'] for row in output}, {'A', 'B'})

        # Only A and B survive: Nc = 151+21, Ne = 11+151.
        pairs = {'A': (151, 11), 'B': (21, 151)}
        expected = {name: proportions_ztest(pair, [172, 162])
                    for name, pair in pairs.items()}
        bh = multipletests([expected[name][1] for name in ('A', 'B')],
                          method='fdr_bh')[1]

        for row in output:
            name = row['gene_name']
            z, p = expected[name]
            self.assertEqual(int(row['gene_length']), 1000)
            self.assertEqual(int(row['N_contacts_plus_1']), pairs[name][0])
            self.assertEqual(int(row['N_counts_RNAseq_plus_1']), pairs[name][1])
            self.assertAlmostEqual(float(row['chP']), z, places=12)
            self.assertAlmostEqual(float(row['pval']), p, places=12)
            self.assertAlmostEqual(float(row['fdr_bh']), bh[('A', 'B').index(name)], places=12)
        self.assertEqual([row['gene_name'] for row in output],
                         sorted(['A', 'B'], key=lambda name: (bh[('A', 'B').index(name)], name)))

    def test_empty_result_writes_header_only_and_no_png(self):
        """No eligible genes must not leave an outdated image."""
        self.out.mkdir(parents=True)
        image = self.out / 'chP.png'
        image.write_bytes(b'outdated plot')
        result = self.run_script('counts', '--N_contacts_min', '1000000')
        self.assertEqual(result.returncode, 0, result.stderr)
        columns, rows = read_tsv(self.out / 'chP.tsv')
        self.assertEqual(columns, OUTPUT_COLUMNS)
        self.assertEqual(rows, [])
        self.assertFalse(image.exists())

    def test_one_eligible_gene_is_safe(self):
        """A single eligible gene has no defined two-proportion comparison."""
        result = self.run_script('counts', '--N_contacts_min', '30')
        self.assertEqual(result.returncode, 0, result.stderr)
        _, rows = read_tsv(self.out / 'chP.tsv')
        self.assertEqual(rows, [])

    def test_strict_contact_threshold(self):
        """The counts cutoff is strict and applied before the +1 adjustment."""
        result = self.run_script('counts', '--N_contacts_min', '20')
        self.assertEqual(result.returncode, 0, result.stderr)
        _, rows = read_tsv(self.out / 'chP.tsv')
        self.assertEqual(rows, [])  # A alone survives, no Z statistic

    def test_reject_invalid_fdr(self):
        """Invalid FDR thresholds should fail at argument parsing."""
        for value in ('0', '-0.1', '1.1'):
            with self.subTest(value=value):
                result = self.run_script('counts', '--fdr_threshold', value)
                self.assertNotEqual(result.returncode, 0)

    def test_png_creation_and_signature(self):
        """Both biotype groups should generate a real PNG."""
        # Multiple significant points per biotype also exercise marginal KDEs.
        write_tsv(self.rna, HEADER, [
            ['A', 'chr1', 1, 1000, 'protein_coding', 'test', 10],
            ['B', 'chr1', 1, 1000, 'lncRNA', 'test', 150],
            ['E', 'chr1', 1, 1000, 'protein_coding', 'test', 300],
            ['F', 'chr1', 1, 1000, 'lncRNA', 'test', 5],
        ])
        write_tsv(self.contacts, HEADER, [
            ['A', 'chr1', 1, 1000, 'protein_coding', 'test', 150],
            ['B', 'chr1', 1, 1000, 'lncRNA', 'test', 20],
            ['E', 'chr1', 1, 1000, 'protein_coding', 'test', 30],
            ['F', 'chr1', 1, 1000, 'lncRNA', 'test', 170],
        ])
        result = self.run_script(plot=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        image = self.out / 'chP.png'
        self.assertTrue(image.is_file())
        self.assertGreater(image.stat().st_size, 5000)
        self.assertEqual(image.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')

    def test_unique_dna_mode_and_filename(self):
        """UU counting should feed both contact-count and chP outputs."""
        write_tsv(self.uu, CONTACT_COLUMNS,
                  [contact(f'a{i}', 'A', 'protein_coding', 100000) for i in range(5)] +
                  [contact(f'b{i}', 'B', 'lncRNA', 100000) for i in range(3)])
        result = self.run_script('unique_dna', '--contacts-uu', self.uu,
                                 '--dist', '500')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.out / 'contact_counts_unique_dna_dist_500.tsv').exists())
        columns, rows = read_tsv(self.out / 'chP_unique_dna_dist_500.tsv')
        self.assertEqual(columns, OUTPUT_COLUMNS)
        self.assertEqual({row['gene_name'] for row in rows}, {'A', 'B'})
        self.assertTrue((self.out / 'chP_unique_dna_dist_500.png').exists())


if __name__ == '__main__':
    unittest.main()
