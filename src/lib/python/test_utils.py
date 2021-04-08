# coding: utf-8

import os
import unittest
import shutil
import tempfile

class TestUtilMethods(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate_sample_list(self):
        from .utils import generate_sample_list_from_samplesheet

        data = generate_sample_list_from_samplesheet("tests/samplesheets/samplesheet1.csv")
        self.assertEqual(data, {
            "29948": {"counter": 1, "projectpath": ""},
            "D21-00844": {"counter": 2, "projectpath": ""},
            "29949": {"counter": 3, "projectpath": ""},
            "D21-00845": {"counter": 4, "projectpath": ""}})

        data = generate_sample_list_from_samplesheet("tests/samplesheets/samplesheet2.csv")
        self.assertEqual(data, {
            "29948": {"counter": 1, "projectpath": "WP3/"},
            "D21-00844": {"counter": 2, "projectpath": "WP2/"},
            "29949": {"counter": 3, "projectpath": "WP3/"},
            "D21-00845": {"counter": 4, "projectpath": "WP2/"}})
