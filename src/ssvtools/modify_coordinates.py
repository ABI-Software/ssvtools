"""
Utility functions for modifying coordinates of SPARC subject-specific nerve scaffolds,
including adopting template trunk coordinates.
"""
import logging
import math
import os
import re

from cmlibs.maths.vectorops import normalize, cross, add, mult, magnitude, set_magnitude
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.fieldmodule import Fieldmodule
from ssvtools.query_structure import get_vagus_trunk_group


logger = logging.getLogger(__name__)


def adopt_template_trunk_coordinates(region, coordinate_field_name, template_region, template_coordinate_field_name,
                                     trunk_group_name, unit_conversion_factor):
    """
    Adopt trunk coordinates from a template vagus scaffold (for example, one derived from the 3D whole body)
    for a subject-specific vagus (SSV) scaffold. For this to work, the SSV vagus trunk needs to have the same number of
    elements along the trunk of a template vagus scaffold. The coordinates field of the template trunk is written to a
    memory buffer in a template region and loaded back into the region to replace the field of the SSV trunk with the
    same name as the template coordinate field. The values are then used to replace the values in the choosen SSV
    coordinate field. This means that the chosen coordinate field of the SSV trunk now follows the path of the template
    trunk. The radius along the resulting trunk is then scaled to the radius along the original SSV trunk using a unit
    conversion factor. As the SSV now has a rescaled trunk, the new position of the SSV branches can be determined from
    the new trunk coordinates.
    :param region: region where SSV is loaded.
    :param coordinate_field_name: name of the chosen coordinate field for SSV.
    :param template_region: template region where the template vagus scaffold is loaded.
    :param template_coordinate_field_name: name of the coordinate field for template vagus scaffold.
    :param trunk_group_name: name of trunk group.
    :param unit_conversion_factor: Factor to bring SSV into the same scale as template scaffold. If none, mesh integral
    of the respective trunk lengths will be used to calculate the required scale factor.
    """
    fieldmodule = region.getFieldmodule()
    mesh3d = fieldmodule.findMeshByDimension(3)
    fieldcache = fieldmodule.createFieldcache()
    coordinates_field = fieldmodule.findFieldByName(coordinate_field_name).castFiniteElement()
    tmp_coordinates_field = fieldmodule.findFieldByName(template_coordinate_field_name).castFiniteElement()

    vagus_trunk_group = get_vagus_trunk_group(fieldmodule)
    trunk_mesh_group3d = vagus_trunk_group.getMeshGroup(mesh3d)
    elements_along_data_trunk = trunk_mesh_group3d.getSize()

    template_fieldmodule = template_region.getFieldmodule()
    template_mesh3d = template_fieldmodule.findMeshByDimension(3)
    vagus_trunk_group = get_vagus_trunk_group(template_fieldmodule)
    trunk_mesh_group3d = vagus_trunk_group.getMeshGroup(template_mesh3d)
    elements_along_template_trunk = trunk_mesh_group3d.getSize()
    assert elements_along_template_trunk == elements_along_data_trunk, \
        "adopt_template_trunk_coordinates. Number of trunk elements in data does not match number of trunk elements in " \
        "template vagus."

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

    # Read in coordinates of the nodes in template trunk group to the region
    sir = template_region.createStreaminformationRegion()
    srm = sir.createStreamresourceMemory()
    sir.setResourceDomainTypes(srm, Field.DOMAIN_TYPE_NODES)
    sir.setResourceGroupName(srm, trunk_group_name)
    sir.setResourceFieldNames(srm, template_coordinate_field_name)
    if unit_conversion_factor is None:
        template_coordinates_field = template_fieldmodule.findFieldByName(
            template_coordinate_field_name).castFiniteElement()
        unit_conversion_factor = getUnitConversionFactor(trunk_group_name, fieldmodule, coordinates_field,
                                                         template_fieldmodule, template_coordinates_field)
    template_region.write(sir)
    result, buffer = srm.getBuffer()
    sir = region.createStreaminformationRegion()
    srm = sir.createStreamresourceMemoryBuffer(buffer)
    result = region.read(sir)

    del srm
    del sir
    del template_mesh3d
    del template_fieldmodule
    del template_region

    derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
    derivative_xi3 = mesh3d.getChartDifferentialoperator(1, 3)

    elem_iter = mesh3d.createElementiterator()
    element = elem_iter.next()
    element_id = element.getIdentifier()
    count = 0
    while element.isValid():
        eft = element.getElementfieldtemplate(coordinates_field, -1)
        # Loop through coordinates field of trunk elements and set values to chosen SSV coordinates field. Also scale
        # radius on template trunk to match relative radius in SSV trunk
        if element_id <= elements_along_template_trunk:
            if element_id == 1:
                ln = [1, 2]
            else:
                ln = [2]
            for i in range(len(ln)):
                node = element.getNode(eft, ln[i])
                fieldcache.setNode(node)
                # get radius and rate of change of radius of template trunk
                _, bx = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, bd1 = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                _, bd2 = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                _, bd3 = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                _, bd12 = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
                _, bd13 = tmp_coordinates_field.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, bx)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, bd1)
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                    set_magnitude(bd2, unit_conversion_factor * bd2_mag_list[count]))
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                    set_magnitude(bd3, unit_conversion_factor * bd3_mag_list[count]))
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1,
                                                    set_magnitude(bd12, unit_conversion_factor * bd12_mag_list[count]))
                coordinates_field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1,
                                                    set_magnitude(bd13, unit_conversion_factor * bd13_mag_list[count]))
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
        element_id = element.getIdentifier()


def getUnitConversionFactor(trunk_group_name, fieldmodule, coordinates_field, template_fieldmodule,
                            template_coordinates_field):
    """
    Derive a unit conversion factor by calculating the scale difference between the length of the SSV scaffold and the
    template trunk length.
    """
    trunk_group = fieldmodule.findFieldByName(trunk_group_name).castGroup()
    centroid_group = fieldmodule.findFieldByName("vagus centroid").castGroup()
    mesh = create_intersection_mesh_group(trunk_group, centroid_group, dimension=1)
    SSV_trunk_length = get_mesh_integral_one(mesh, coordinates_field)

    trunk_group = template_fieldmodule.findFieldByName(trunk_group_name).castGroup()
    centroid_group = template_fieldmodule.findFieldByName("vagus centroid").castGroup()
    mesh = create_intersection_mesh_group(trunk_group, centroid_group, dimension=1)
    template_trunk_length = get_mesh_integral_one(mesh, template_coordinates_field)

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

def get_mesh_integral_one(mesh, coordinates, number_of_points=4):
    """
    Calculate the length of a mesh using meshIntegral.
    """
    fieldmodule = mesh.getFieldmodule()
    with ChangeManager(fieldmodule):
        mesh_integral = \
            fieldmodule.createFieldMeshIntegral(fieldmodule.createFieldConstant(1.0), coordinates, mesh)
        mesh_integral.setNumbersOfPoints([number_of_points])
        fieldcache = fieldmodule.createFieldcache()
        result, value = mesh_integral.evaluateReal(fieldcache, 1)
        del mesh_integral
    return value
