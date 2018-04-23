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
import hicap


class CommandExecuteTestCase(unittest.TestCase):

    def test_execute_command_1(self):
        result = hicap.execute_command('echo -n test')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test')
        self.assertEqual(result.stderr, '')

    def test_execute_command_2(self):
        result = hicap.execute_command('invalid_command', check=False)
        self.assertEqual(result.returncode, 127)
        self.assertEqual(result.stdout, '')

    def test_execute_command_3(self):
        with self.assertRaises(SystemExit):
            hicap.execute_command('invalid_command')


class FastaParserTestCase(unittest.TestCase):

    def test_read_query_fasta_1(self):
        input_fp = pathlib.Path(tests_directory, 'data/good.fasta')
        fasta = hicap.read_query_fasta(input_fp)
        self.assertEqual(fasta['good'], 'atgc')

    def test_read_query_fasta_2(self):
        input_fp = pathlib.Path(tests_directory, 'data/bad.fasta')
        with self.assertRaises(SystemExit):
            hicap.read_query_fasta(input_fp)
