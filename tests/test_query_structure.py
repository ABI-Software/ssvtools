from cmlibs.zinc.context import Context
from cmlibs.zinc.result import RESULT_OK
from ssvtools.modify_coordinates import adopt_template_trunk_coordinates
from ssvtools.query_structure import evaluate_branch_start_coordinates, get_marker_data, get_vagus_structure_maps, \
    get_vagus_trunk_group
from testutils import assertAlmostEqualList
import os
import unittest


here = os.path.abspath(os.path.dirname(__file__))


class SSVToolsTestCase(unittest.TestCase):

    def test_query_branches(self):
        """
        Get positions and orientations of branches from a test subject-specific-vagus SSV.
        Uses a simple test scaffold, but any vagus scaffold output by SPARC Mapping Tools
        Scaffold Creator, including from REVA data should work.
        """
        data_file = os.path.join(here, "resources", "vagus_test_scaffold1.exf")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertEqual(RESULT_OK, region.readFile(data_file))
        fieldmodule = region.getFieldmodule()
        structure_map, common_branch_map = get_vagus_structure_maps(fieldmodule)
        expected_structure_map = {
            'left vagus nerve': [None, ['left A thoracic cardiopulmonary branch of vagus nerve',
                                        'left B thoracic cardiopulmonary branch of vagus nerve',
                                        'left superior laryngeal nerve']],
            'left A branch of superior laryngeal nerve': ['left superior laryngeal nerve', []],
            'left A thoracic cardiopulmonary branch of vagus nerve': ['left vagus nerve', []],
            'left B thoracic cardiopulmonary branch of vagus nerve': ['left vagus nerve', []],
            'left superior laryngeal nerve': ['left vagus nerve', ['left A branch of superior laryngeal nerve']]}
        self.assertEqual(expected_structure_map, structure_map)
        expected_common_branch_map = {
            'left branch of superior laryngeal nerve': [
                'left A branch of superior laryngeal nerve'],
            'left thoracic cardiopulmonary branch of vagus nerve': [
                'left A thoracic cardiopulmonary branch of vagus nerve',
                'left B thoracic cardiopulmonary branch of vagus nerve']}
        self.assertEqual(common_branch_map, expected_common_branch_map)
        vagus_coordinates = fieldmodule.findFieldByName('vagus coordinates')
        branch_start_coordinates = evaluate_branch_start_coordinates(
            vagus_coordinates, structure_map['left vagus nerve'][1])
        expected_branch_start_coordinates = [
            ('left A thoracic cardiopulmonary branch of vagus nerve',
                [-0.0003328408968624165, 1.7493737182796838e-05, 0.3813658363063589],
                [-0.9106917613635638, 0.38295521752225525, 0.15487355215469165]),
            ('left B thoracic cardiopulmonary branch of vagus nerve',
                [0.0010020650864467912, -0.001554513928071552, 0.40674113258208205],
                [0.5645932995214857, -0.630854592172089, 0.5322188362608271]),
            ('left superior laryngeal nerve',
                [0.0008767095326078595, 0.00045882124897275354, 0.1315331326187551],
                [0.6662449900204912, 0.64287539363479, -0.3779270320200837])]
        self.assertEqual(len(branch_start_coordinates), len(expected_branch_start_coordinates))
        TOL = 1.0E-8
        for i in range(len(branch_start_coordinates)):
            expected_branch_name, expected_start_location, expected_direction = expected_branch_start_coordinates[i]
            branch_name, start_location, direction = branch_start_coordinates[i]
            self.assertEqual(branch_name, expected_branch_name)
            assertAlmostEqualList(self, start_location, expected_start_location, TOL)
            assertAlmostEqualList(self, direction, expected_direction, TOL)
        straight_coordinates = fieldmodule.findFieldByName('straight coordinates')
        marker_data = get_marker_data(fieldmodule, straight_coordinates)
        expected_marker_data = [
            ('left level of superior border of jugular foramen on the vagus nerve',
             (1, [0.6843240000000003, 0.5, 0.5]),
             [0.0, 0.0, 0.02737296], [-6.50506937191939e-14, -1.310900721441712e-14, 1774.9571795976287]),
            ('left level of inferior border of jugular foramen on the vagus nerve', (2, [0.1087380000000001, 0.5, 0.5]),
             [0.0, 0.0, 0.04434952], [1.0209665570151899e-13, 9.768198918855594e-14, 2991.9103185208987]),
            ('left level of angle of the mandible on the vagus nerve', (4, [0.1332685000000003, 0.5, 0.5]),
             [0.0, 0.0, 0.12533074], [1.4483477432177637e-14, 2.1376726263214913e-14, 8245.63333162319]),
            ('left level of carotid bifurcation on the vagus nerve', (4, [0.9345910000000012, 0.5, 0.5]),
             [0.0, 0.0, 0.15738364], [-8.404728289972044e-14, 1.2993009323710421e-14, 10174.73968280258]),
            ('left level of laryngeal prominence on the vagus nerve', (6, [0.1354834999999998, 0.5, 0.5]),
             [0.0, 0.0, 0.20541934], [3.3090283510828584e-14, -1.6217927879715146e-14, 13115.836142905777]),
            ('left level of superior border of the clavicle on the vagus nerve', (9, [0.4619939999999999, 0.5, 0.5]),
             [0.0, 0.0, 0.33847976], [3.197424311348175e-14, 3.8604867982647756e-14, 21182.365707090583]),
            ('left level of jugular notch on the vagus nerve', (10, [0.515577750000001, 0.5, 0.5]), [0.0, 0.0, 0.38062311],
             [4.9203238902202615e-14, 5.4205666260696657e-14, 23763.36090724233]),
            ('left level of sternal angle on the vagus nerve', (13, [0.09881600000000146, 0.5, 0.5]),
             [0.0, 0.0, 0.48395264], [7.528280380809436e-14, 7.334399515678671e-14, 30050.97537892891])]
        STOL = 1.0E-4
        for i in range(len(expected_marker_data)):
            expected_marker_name, expected_marker_location, expected_marker_vagus_coordinates, \
                expected_marker_straight_coordinates = expected_marker_data[i]
            marker_name, marker_location, marker_vagus_coordinates, marker_straight_coordinates = marker_data[i]
            self.assertEqual(marker_name, expected_marker_name)
            self.assertEqual(marker_location[0], expected_marker_location[0])
            assertAlmostEqualList(self, marker_location[1], expected_marker_location[1], TOL)
            assertAlmostEqualList(self, marker_vagus_coordinates, expected_marker_vagus_coordinates, TOL)
            assertAlmostEqualList(self, marker_straight_coordinates, expected_marker_straight_coordinates, STOL)

    def test_straight_to_geometric(self):
        """
        Modify the straight coordinates from a test subject-specific-vagus (SSV) by adopting its geometric coordinates
        along the trunk. Uses a simple test scaffold, but any vagus scaffold output by SPARC Mapping Tools Scaffold
        Creator, including from REVA data should work.
        """
        data_file_path = os.path.join(here, "resources", "vagus_test_scaffold2.exf")
        unit_conversion_factor = 1.0
        trunk_group_name = "left vagus nerve"
        context = Context("Test")
        region = context.getDefaultRegion()
        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, region.readFile(data_file_path))
        coordinates_field_name = "straight coordinates"
        coordinates_field = fieldmodule.findFieldByName(coordinates_field_name).castFiniteElement()

        template_region = region.createRegion()
        template_fieldmodule = template_region.getFieldmodule()
        self.assertEqual(RESULT_OK, template_region.readFile(data_file_path))
        template_coordinates_field_name = "coordinates"
        template_coordinates_field = template_fieldmodule.findFieldByName(template_coordinates_field_name).castFiniteElement()

        adopt_template_trunk_coordinates(region, coordinates_field, template_region, template_coordinates_field,
                                         trunk_group_name, unit_conversion_factor)

        # region.writeFile('straight_to_geometric.exf')
        coordinates_field = fieldmodule.findFieldByName(coordinates_field_name).castFiniteElement()
        structure_map, common_branch_map = get_vagus_structure_maps(fieldmodule)
        vagus_trunk_group = get_vagus_trunk_group(fieldmodule)
        vagus_trunk_group_name = vagus_trunk_group.getName()
        trunk_branch_list = structure_map[vagus_trunk_group_name][1]
        branch_start_data = evaluate_branch_start_coordinates(coordinates_field, trunk_branch_list)

        expected_branch_start_data = \
            [('left A thoracic cardiopulmonary branch of vagus nerve',
              [20634.120654634753, -2944.00512888507, -607.5696487203066],
              [0.12153923672191891, -0.9919656891879115, -0.035103923014460706]),
             ('left B thoracic cardiopulmonary branch of vagus nerve',
              [22162.25852383412, -3220.7718823220207, -620.3932944321155],
              [0.7452977441582963, 0.6546665028595009, -0.12626576173399476]),
             ('left superior laryngeal nerve', [5920.133012880977, -4451.293552269674, -195.3146602807042],
              [-0.8472154614382229, 0.5304071806950975, -0.029902919039740395])]

        TOL = 1.0E-8
        for i in range(len(expected_branch_start_data)):
            expected_branch_name, expected_start_x, expected_direction = expected_branch_start_data[i]
            branch_name, start_x, direction = branch_start_data[i]
            self.assertEqual(branch_name, expected_branch_name)
            assertAlmostEqualList(self, start_x, expected_start_x, TOL)
            assertAlmostEqualList(self, direction, expected_direction, TOL)

    def test_geometric_to_straight(self):
        """
        Modify the geometric coordinates from a test subject-specific-vagus (SSV) by adopting its straight coordinates
        along the trunk. Also, used None for unit_conversion_factor in this test case. Uses a simple test scaffold, but
        any vagus scaffold output by SPARC Mapping Tools Scaffold Creator, including from REVA data should work.
        """
        data_file_path = os.path.join(here, "resources", "vagus_test_scaffold2.exf")
        unit_conversion_factor = None
        trunk_group_name = "left vagus nerve"
        context = Context("Test")
        region = context.getDefaultRegion()
        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, region.readFile(data_file_path))
        coordinates_field_name = "coordinates"
        coordinates_field = fieldmodule.findFieldByName(coordinates_field_name).castFiniteElement()

        template_region = region.createRegion()
        template_fieldmodule = template_region.getFieldmodule()
        self.assertEqual(RESULT_OK, template_region.readFile(data_file_path))
        template_coordinates_field_name = "straight coordinates"
        template_coordinates_field = \
            template_fieldmodule.findFieldByName(template_coordinates_field_name).castFiniteElement()

        adopt_template_trunk_coordinates(region, coordinates_field, template_region, template_coordinates_field,
                                         trunk_group_name, unit_conversion_factor)

        # region.writeFile('geometric_to_straight.exf')
        coordinates_field = fieldmodule.findFieldByName(coordinates_field_name).castFiniteElement()
        structure_map, common_branch_map = get_vagus_structure_maps(fieldmodule)
        vagus_trunk_group = get_vagus_trunk_group(fieldmodule)
        vagus_trunk_group_name = vagus_trunk_group.getName()
        trunk_branch_list = structure_map[vagus_trunk_group_name][1]
        branch_start_data = evaluate_branch_start_coordinates(coordinates_field, trunk_branch_list)

        expected_branch_start_data = \
            [('left A thoracic cardiopulmonary branch of vagus nerve',
              [-11.396312835632429, 0.046011947444071666, 23760.81526778722],
              [-0.8902562546591324, 0.36093917883438575, 0.27778896706430495]), (
             'left B thoracic cardiopulmonary branch of vagus nerve',
             [29.374459368699924, -51.59767903274964, 25311.09739308472],
             [0.5373745137850673, -0.6508445619590075, 0.536311465570754]), (
             'left superior laryngeal nerve', [46.79098252744873, 18.715786185700324, 8593.065332472212],
             [0.7403494295872118, 0.6408371734692249, -0.20300354482093294])]

        TOL = 1.0E-8
        for i in range(len(expected_branch_start_data)):
            expected_branch_name, expected_start_x, expected_direction = expected_branch_start_data[i]
            branch_name, start_x, direction = branch_start_data[i]
            self.assertEqual(branch_name, expected_branch_name)
            assertAlmostEqualList(self, start_x, expected_start_x, TOL)
            assertAlmostEqualList(self, direction, expected_direction, TOL)


if __name__ == "__main__":
    unittest.main()
