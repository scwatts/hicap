'''
Copyright 2018 Stephen Watts
https://github.com/scwatts/hicap

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


import logging
import pathlib
import unittest


from . import tests_directory
import hicap.blast


class CreateDatabaseTestCase(unittest.TestCase):

    def setUp(self):
        self.blast_db_fp = None
        self.blast_db_ext = ('nhr', 'nin', 'nsq')

    def tearDown(self):
        for ext in self.blast_db_ext:
            blast_db_fp_part = pathlib.Path('%s.%s' % (self.blast_db_fp, ext))
            if blast_db_fp_part.exists():
                blast_db_fp_part.unlink()

    def test_create_database_1(self):
        fasta_fp = pathlib.Path(tests_directory, 'data/good.fasta')
        self.blast_db_fp = hicap.blast.create_blast_database(fasta_fp, tests_directory)
        for ext in self.blast_db_ext:
            blast_db_fp_part = pathlib.Path('%s.%s' % (self.blast_db_fp, ext))
            self.assertTrue(blast_db_fp_part.exists())

    def test_create_database_2(self):
        fasta_fp = pathlib.Path(tests_directory, 'data/bad.fasta')
        with self.assertRaises(SystemExit):
            hicap.blast.create_blast_database(fasta_fp, tests_directory)


class QueryDatabaseTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        database_fp = pathlib.Path(tests_directory, 'data/type_a.fasta')
        cls.blast_db_fp = hicap.blast.create_blast_database(database_fp, tests_directory)
        cls.blast_db_ext = ('nhr', 'nin', 'nsq')
        hicap.blast.BlastFormat = {'qseqid': str,
                                   'sseqid': str,
                                   'qlen': int,
                                   'slen': int,
                                   'qstart': int,
                                   'qend': int,
                                   'sstart': int,
                                   'send': int,
                                   'length': int,
                                   'mismatch': int,
                                   'gaps': int}

    @classmethod
    def tearDownClass(cls):
        for ext in cls.blast_db_ext:
            blast_db_fp_part = pathlib.Path('%s.%s' % (cls.blast_db_fp, ext))
            blast_db_fp_part.unlink()

    def test_query_database_1(self):
        first_hit = 'Z37516\tZ37516\t5882\t5882\t2161\t5882\t2161\t5882\t3722\t0\t0'
        last_hit = 'Z37516\tZ37516\t5882\t5882\t2965\t2976\t5696\t5685\t12\t0\t0'
        query_fp = pathlib.Path(tests_directory, 'data/type_a_variant.fasta')
        blast_stdout = hicap.blast.align(query_fp, self.blast_db_fp)
        hits = [l for l in blast_stdout.rstrip().split('\n')]
        self.assertEqual(hits[0], first_hit)
        self.assertEqual(hits[-1], last_hit)

    def test_query_database_2(self):
        query_fp = pathlib.Path(tests_directory, 'data/good.fasta')
        blast_stdout = hicap.blast.align(query_fp, self.blast_db_fp)
        self.assertEqual(blast_stdout, '')
