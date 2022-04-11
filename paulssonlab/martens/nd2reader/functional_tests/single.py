"""
These tests require that you have a specific ND2 file created by the developer of nd2reader. You will never need to
run them unless you're Jim Rybarski.

"""
from nd2reader import Nd2
from datetime import datetime
import unittest


class SingleTests(unittest.TestCase):
    def setUp(self):
        self.nd2 = Nd2("/var/nd2s/single.nd2")

    def tearDown(self):
        self.nd2.close()

    def test_shape(self):
        self.assertEqual(self.nd2.height, 512)
        self.assertEqual(self.nd2.width, 512)

    def test_date(self):
        self.assertEqual(self.nd2.date, datetime(2015, 10, 15, 9, 33, 5))

    def test_length(self):
        self.assertEqual(len(self.nd2), 1)

    def test_channels(self):
        self.assertEqual(self.nd2.channels, ["Quad Band 2"])

    def test_pixel_size(self):
        self.assertGreater(self.nd2.pixel_microns, 0.0)

    def test_actual_length(self):
        count = 0
        for image in self.nd2:
            if image is not None:
                count += 1
        self.assertEqual(len(self.nd2), count)

    def test_frames(self):
        self.assertEqual(len(self.nd2.frames), 1)

    def test_fovs(self):
        self.assertEqual(len(self.nd2.fields_of_view), 1)

    def test_z_levels(self):
        self.assertTupleEqual(tuple(self.nd2.z_levels), (0,))

    def test_image(self):
        image = self.nd2[0]
        self.assertIsNotNone(image)

    def test_iteration(self):
        images = [image for image in self.nd2]
        self.assertEqual(len(images), 1)

    def test_iteration_step(self):
        images = [image for image in self.nd2[::2]]
        self.assertEqual(len(images), 1)

    def test_iteration_backwards(self):
        images = [image for image in self.nd2[::-1]]
        self.assertEqual(len(images), 1)

    def test_select_bounds_wrong(self):
        images = [i for i in self.nd2.select(start=0, stop=12481247)]
        self.assertEqual(len(images), 1)
