"""
Utility functions for modifying coordinates of SPARC subject-specific nerve scaffolds,
including adopting template trunk coordinates.
"""
import logging
import math

from cmlibs.maths.vectorops import normalize, cross, add, mult, magnitude, set_magnitude
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.fieldmodule import Fieldmodule
from ssvtools.query_structure import get_vagus_trunk_group


logger = logging.getLogger(__name__)


def adopt_template_trunk_coordinates(region, coordinates_field, template_region, template_coordinates_field,
                                     trunk_group_name, unit_conversion_factor):
    """
    Adopt trunk coordinates from a template vagus scaffold (for example, one derived from the 3D whole body)
    for a subject-specific vagus (SSV) scaffold. For this to work, the SSV vagus trunk needs to have the same number of
    elements along the trunk of a template vagus scaffold. The chosen coordinate field of the SSV trunk now follows the
    path of the template trunk. For radius, the trunk uses the radius of the SSV but converts units using the
    unit_conversion_factor. For branches, the result is SSV branches in their equivalent directions but now relative
    to the template trunk coordinates; these radiate in straight lines from the new trunk, as do branches of branches.
    Note that SSV and template vagus nerve should be for the same side, left or right; a logger warning is output
    if this is not the case, but otherwise the function proceeds.
    :param region: region where SSV is loaded.
    :param coordinates_field: Chosen coordinates field for SSV. Must be a field from the SSV region.
    :param template_region: template region where the template vagus scaffold is loaded.
    :param template_coordinates_field: Coordinates field for template vagus scaffold. Must be a field from the
    template_region.
    :param trunk_group_name: name of trunk group.
    :param unit_conversion_factor: Factor to bring SSV into the same scale as template scaffold. If None, the ratio of
    trunk lengths is used to calculate a unit conversion factor as a power of 10.
    """
    fieldmodule = region.getFieldmodule()
    mesh3d = fieldmodule.findMeshByDimension(3)
    fieldcache = fieldmodule.createFieldcache()

    vagus_trunk_group = get_vagus_trunk_group(fieldmodule)
    trunk_mesh_group3d = vagus_trunk_group.getMeshGroup(mesh3d)
    trunk_group_name = vagus_trunk_group.getName()
    elements_along_data_trunk = trunk_mesh_group3d.getSize()

    template_fieldmodule = template_region.getFieldmodule()
    template_mesh3d = template_fieldmodule.findMeshByDimension(3)
    template_fieldcache = template_fieldmodule.createFieldcache()

    template_vagus_trunk_group = get_vagus_trunk_group(template_fieldmodule)
    template_trunk_group_name = template_vagus_trunk_group.getName()
    if trunk_group_name != template_trunk_group_name:
        logger.warning("adopt_template_trunk_coordinates. SSV and template have different trunk group names: " +
                       trunk_group_name + " vs. " + template_trunk_group_name)

    if ('left' in trunk_group_name and 'right' in template_trunk_group_name) or \
            ('right' in trunk_group_name and 'left' in template_trunk_group_name) :
        logger.warning("adopt_template_trunk_coordinates. SSV and template have trunks from different side: " +
                       trunk_group_name + " vs. " + template_trunk_group_name)

    template_trunk_mesh_group3d = template_vagus_trunk_group.getMeshGroup(template_mesh3d)
    elements_along_template_trunk = template_trunk_mesh_group3d.getSize()
    assert elements_along_template_trunk == elements_along_data_trunk, \
        "adopt_template_trunk_coordinates. Number of trunk elements in data does not match number of trunk elements " \
        "in template vagus."

    # Get radius from data trunk
    elem_iter = mesh3d.createElementiterator()
    element = elem_iter.next()
    element_id = element.getIdentifier()
    bd2_mag_list = []
    bd3_mag_list = []
    bd12_mag_list = []
    bd13_mag_list = []
    ln = 2

    # Loop through SSV trunk elements to store radius information
    while element.isValid() and element_id <= elements_along_template_trunk:
        eft = element.getElementfieldtemplate(coordinates_field, -1)
        if element_id == 1:
            ln = [1, 2]
        else:
            ln = [2]
        for i in range(len(ln)):
            node = element.getNode(eft, ln[i])
            fieldcache.setNode(node)
            _, bd2 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            _, bd3 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            _, bd12 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
            _, bd13 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
            bd2_mag_list.append(magnitude(bd2))
            bd3_mag_list.append(magnitude(bd3))
            bd12_mag_list.append(magnitude(bd12))
            bd13_mag_list.append(magnitude(bd13))
        element = elem_iter.next()
        element_id = element.getIdentifier()

    # Read in coordinates of the nodes in template trunk group to the template region
    sir = template_region.createStreaminformationRegion()
    srm = sir.createStreamresourceMemory()
    sir.setResourceDomainTypes(srm, Field.DOMAIN_TYPE_NODES)
    sir.setResourceGroupName(srm, trunk_group_name)
    template_coordinate_field_name = template_coordinates_field.getName()
    sir.setResourceFieldNames(srm, template_coordinate_field_name)
    if unit_conversion_factor is None:
        template_coordinates_field = template_fieldmodule.findFieldByName(
            template_coordinate_field_name).castFiniteElement()
        unit_conversion_factor = get_unit_conversion_factor(trunk_group_name, fieldmodule, coordinates_field,
                                                            template_fieldmodule, template_coordinates_field)
    template_region.write(sir)

    derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
    derivative_xi3 = mesh3d.getChartDifferentialoperator(1, 3)

    # Loop through coordinates field of trunk elements and set values to chosen SSV coordinates field. Also scale
    # radius on template trunk to match relative radius in SSV trunk
    elem_iter = mesh3d.createElementiterator()
    element = elem_iter.next()
    count = 0
    with ChangeManager(fieldmodule):
        while element.isValid():
            element_id = element.getIdentifier()
            template_element = template_mesh3d.findElementByIdentifier(element_id)
            template_eft = template_element.getElementfieldtemplate(template_coordinates_field, -1)
            eft = element.getElementfieldtemplate(coordinates_field, -1)

            if element_id <= elements_along_template_trunk:
                if element_id == 1:
                    ln = [1, 2]
                else:
                    ln = [2]
                for i in range(len(ln)):
                    # get radius and rate of change of radius of template trunk
                    template_node = template_element.getNode(template_eft, ln[i])
                    template_fieldcache.setNode(template_node)
                    _, bx = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                         Node.VALUE_LABEL_VALUE, 1, 3)
                    _, bd1 = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                          Node.VALUE_LABEL_D_DS1, 1, 3)
                    _, bd2 = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                          Node.VALUE_LABEL_D_DS2, 1, 3)
                    _, bd3 = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                          Node.VALUE_LABEL_D_DS3, 1, 3)
                    _, bd12 = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                           Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
                    _, bd13 = template_coordinates_field.getNodeParameters(template_fieldcache, -1,
                                                                           Node.VALUE_LABEL_D2_DS1DS3, 1, 3)

                    node = element.getNode(eft, ln[i])
                    fieldcache.setNode(node)
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, bx)
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, bd1)
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                        set_magnitude(bd2, unit_conversion_factor *
                                                                      bd2_mag_list[count]))
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                        set_magnitude(bd3, unit_conversion_factor *
                                                                      bd3_mag_list[count]))
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1,
                                                        set_magnitude(bd12, unit_conversion_factor *
                                                                      bd12_mag_list[count]))
                    coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1,
                                                        set_magnitude(bd13, unit_conversion_factor *
                                                                      bd13_mag_list[count]))
                    count += 1
            else:  # Make branches radiating from template trunk
                local_nodes_count = eft.getNumberOfLocalNodes()
                ln = 2
                if local_nodes_count > 2:
                    fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5])  # as xi1 is along the branch
                    _, ax = coordinates_field.evaluateReal(fieldcache, 3)
                    _, ad1 = coordinates_field.evaluateDerivative(derivative_xi1, fieldcache, 3)
                    _, ad3 = coordinates_field.evaluateDerivative(derivative_xi3, fieldcache, 3)
                    # remove shear on remainder of branch:
                    dir1 = normalize(ad1)
                    dir2 = normalize(cross(ad3, dir1))
                    dir3 = cross(dir1, dir2)
                    ln = 3
                else:
                    ax, ad1 = x, d1

                node = element.getNode(eft, ln)
                fieldcache.setNode(node)
                _, bd2 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                _, bd3 = coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

                x = add(ax, ad1)
                d1 = ad1
                d2 = mult(dir2, magnitude(bd2) * unit_conversion_factor)
                d3 = mult(dir3, magnitude(bd3) * unit_conversion_factor)
                d12 = [0.0, 0.0, 0.0]
                d13 = [0.0, 0.0, 0.0]

                node = element.getNode(eft, ln)
                fieldcache.setNode(node)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d12)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, d13)

            element = elem_iter.next()

    del template_fieldmodule
    del template_region


def get_unit_conversion_factor(trunk_group_name, fieldmodule, coordinates_field, template_fieldmodule,
                               template_coordinates_field):
    """
    Derive a unit conversion factor by calculating the scale difference between the length of the SSV scaffold and the
    template trunk length.
    """
    trunk_group = fieldmodule.findFieldByName(trunk_group_name).castGroup()
    centroid_group = fieldmodule.findFieldByName("vagus centroid").castGroup()
    mesh = create_intersection_mesh_group(trunk_group, centroid_group, dimension=1)
    SSV_trunk_length = evaluate_mesh_integral_one(mesh, coordinates_field)

    trunk_group = template_fieldmodule.findFieldByName(trunk_group_name).castGroup()
    centroid_group = template_fieldmodule.findFieldByName("vagus centroid").castGroup()
    mesh = create_intersection_mesh_group(trunk_group, centroid_group, dimension=1)
    template_trunk_length = evaluate_mesh_integral_one(mesh, template_coordinates_field)

    if SSV_trunk_length <= 0 or template_trunk_length <= 0:
        raise ValueError("Values must be positive")

    log_diff = math.log10(SSV_trunk_length) - math.log10(template_trunk_length)
    scale_power = round(log_diff)
    unit_conversion_factor = 1.0 / (10 ** scale_power)

    return unit_conversion_factor


def create_intersection_mesh_group(group1, group2, dimension):
    """
    Create a mesh group which is an intersection between group1 and group2.
    """
    fieldmodule = group1.getFieldmodule()
    with ChangeManager(fieldmodule):
        group = fieldmodule.createFieldGroup()
        mesh = fieldmodule.findMeshByDimension(dimension)
        mesh_group = group.createMeshGroup(mesh)
        mesh_group.addElementsConditional(fieldmodule.createFieldAnd(group1, group2))
    return mesh_group


def evaluate_mesh_integral_one(mesh, coordinates, number_of_points=4):
    """
    Calculate the length, area or volume of the coordinates over the mesh, using Gaussian quadrature.
    :param mesh: Mesh to integrate over, of any dimension.
    :param coordinates: The coordinates field to integrate over.
    :param number_of_points: Number of Gaussian quadrature points. Default 4 is the current maximum supported.
    :return: Length of the mesh if 1D, area if 2D, volume if 3D.
    """
    fieldmodule = mesh.getFieldmodule()
    with ChangeManager(fieldmodule):
        mesh_integral = fieldmodule.createFieldMeshIntegral(
            fieldmodule.createFieldConstant(1.0), coordinates, mesh)
        mesh_integral.setNumbersOfPoints([number_of_points])
        fieldcache = fieldmodule.createFieldcache()
        result, value = mesh_integral.evaluateReal(fieldcache, 1)
        del mesh_integral
    return value

