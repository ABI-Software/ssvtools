"""
Utility functions for querying structure of SPARC subject-specific nerve scaffolds,
including trunk information, branch positions and orientations, and level marker locations.
"""
import logging
import re

from cmlibs.maths.vectorops import normalize
from cmlibs.utils.zinc.field import get_group_list
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.field import Field
from cmlibs.zinc.fieldmodule import Fieldmodule


logger = logging.getLogger(__name__)


# following keyword is present in names of all annotations marking part of the trunk
TRUNK_KEYWORD = 'vagus nerve'
# exception for branches containing the following text:
NON_TRUNK_KEYWORD = 'of vagus nerve'
# following keywords are part of names of all annotations marking branches (without trunk_keywords)
BRANCH_KEYWORDS = ['branch', 'nerve']


def get_vagus_trunk_group(fieldmodule: Fieldmodule):
    """
    Get the Zinc Group containing the vagus trunk part of the model.
    :param fieldmodule: Fieldmodule of Zinc Region containing subject-specific-vagus nerve scaffold.
    :return: Zinc Group for left/right vagus nerve or None if none.
    """
    for side in ('left', 'right'):
        group_name = side + ' ' + TRUNK_KEYWORD
        # a group contains a set of model objects (elements, nodes), but is also
        # a field returning True at any location inside those parts, False outside.
        group = fieldmodule.findFieldByName(group_name).castGroup()
        if group.isValid():
            return group
    return None


def get_vagus_structure_maps(fieldmodule: Fieldmodule):
    """
    Get the structure of the nerve:
    1. map of all nerve tree trunk/single branch group names to parent trunk/branch group name (or None) and
    list of child branch group names (can be empty).
    2. map of all common branch names to the list of variant single branch names A, B, C etc.
    :param fieldmodule: Fieldmodule of Zinc Region containing subject-specific-vagus nerve scaffold.
    :return: dict group name -> (parent group name, list of child group names),
    dict common group name -> list of variant group names.
    """
    group_list = get_group_list(fieldmodule)
    # vagus scaffolds have common groups for all branches with variant names containing letters A, B, C, etc.
    # first determine these groups as a map from common group name to list of variant group names
    common_group_map = {}
    for group in group_list:
        group_name = group.getName()
        common_name = re.sub(r'\b[A-Z]\b\s?', '', group_name).strip()
        if common_name != group_name:
            variant_group_list = common_group_map.get(common_name)
            if not variant_group_list:
                variant_group_list = []
            variant_group_list.append(group_name)
            common_group_map[common_name] = variant_group_list
    trunk_group = get_vagus_trunk_group(fieldmodule)
    mesh3d = fieldmodule.findMeshByDimension(3)
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
    if not (trunk_group and coordinates.isValid()):
        logger.warning('get_vagus_structure_maps: Missing vagus trunk group')
        return {}, {}
    trunk_group_name = trunk_group.getName()
    group_nodeset_map = {trunk_group_name: trunk_group.getNodesetGroup(nodes)}
    single_branch_groups = []
    # map of group name to [parent group name, list of child group names]
    structure_map = {trunk_group_name: [None, []]}
    for group in group_list:
        group_name = group.getName()
        # exclude common groups
        if common_group_map.get(group_name):
            continue
        # exclude trunk groups
        if (TRUNK_KEYWORD in group_name) and not (NON_TRUNK_KEYWORD in group_name):
            continue
        # exclude non-matches to branch keywords
        if not any(branch_keyword in group_name for branch_keyword in BRANCH_KEYWORDS):
            continue
        # exclude groups with no 3-D elements
        if group.getMeshGroup(mesh3d).getSize() <= 0:
            continue
        single_branch_groups.append(group)
        group_nodeset_map[group_name] = group.getNodesetGroup(nodes)
        structure_map[group_name] = [None, []]
    for group in single_branch_groups:
        group_name = group.getName()
        # the first element in each branch starts with 2 nodes in the parent branch
        # whose parameters are general-linearly-mapped to get starting branch parameters
        # determine parent group from this
        element = group.getMeshGroup(mesh3d).createElementiterator().next()
        eft = element.getElementfieldtemplate(coordinates, -1)
        node = element.getNode(eft, 1)
        for compare_group_name, compare_nodeset_group in group_nodeset_map.items():
            if compare_group_name == group_name:
                continue
            if compare_nodeset_group.containsNode(node):
                parent_group_name = compare_group_name
                structure_map[group_name][0] = parent_group_name
                structure_map[parent_group_name][1].append(group_name)
                break
    return structure_map, common_group_map


def evaluate_branch_start_coordinates(coordinate_field: Field, branch_group_names):
    """
    Evaluate the start coordinates and normalized (unit) direction of the listed branches.
    :param coordinate_field: Zinc coordinate Field to evaluate.
    :param branch_group_names: List of branch group names.
    :return: list of (branch_group_names, start_x, start_direction)
    """
    fieldmodule = coordinate_field.getFieldmodule()
    mesh3d = fieldmodule.findMeshByDimension(3)
    derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
    fieldcache = fieldmodule.createFieldcache()
    branch_start_coordinates = []
    for branch_group_name in branch_group_names:
        group = fieldmodule.findFieldByName(branch_group_name).castGroup()
        mesh_group = group.getMeshGroup(mesh3d)
        element = mesh_group.createElementiterator().next()  # get the first element
        # this is the centre of the branch at its root
        fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5])
        _, start_x = coordinate_field.evaluateReal(fieldcache, 3)
        _, start_d1 = coordinate_field.evaluateDerivative(derivative_xi1, fieldcache, 3)
        direction = normalize(start_d1)
        branch_start_coordinates.append((branch_group_name, start_x, direction))
    return branch_start_coordinates


def get_marker_data(fieldmodule, host_coordinate_field=None):
    """
    Get the list of vagus level marker points down the trunk, in order from top to bottom.
    :param fieldmodule: 
    :param host_coordinate_field: Optional host coordianates to evaluated e.g. straight coordinates
    or coordinates; vagus coordinates are already known to the markers.
    :return: List of (marker name, marker location (element number, [element local coordinates 1, 2, 3]),
                      marker vagus coordinates [, optional host coordinates])
    """
    marker_group = fieldmodule.findFieldByName("marker").castGroup()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    marker_nodeset_group = marker_group.getNodesetGroup(nodes)
    if marker_nodeset_group.getSize() <= 0:
        logger.warning('get_marker_data: No level markers')
        return []
    # Zinc is an interactive system which notifies clients of changes to fields etc.
    # The following caches multiple changes to minimise expensive messaging
    marker_data = []
    with ChangeManager(fieldmodule):
        marker_name = fieldmodule.findFieldByName("marker_name")
        marker_location = fieldmodule.findFieldByName("marker_location")
        # markers store their vagus coordinates in a separate field from the scaffold presently
        marker_vagus_coordinates = fieldmodule.findFieldByName("marker vagus coordinates")
        # a common idiom in Zinc is to use its rich field operators to make fields returning the values we want
        # here we evaluate the host field at the location stored at the marker
        marker_host_coordinates = fieldmodule.createFieldEmbedded(
            host_coordinate_field, marker_location) if host_coordinate_field else None
        fieldcache = fieldmodule.createFieldcache()
        nodeiterator = marker_nodeset_group.createNodeiterator()
        node = nodeiterator.next()
        while node.isValid():
            fieldcache.setNode(node)
            name = marker_name.evaluateString(fieldcache)
            element, xi = marker_location.evaluateMeshLocation(fieldcache, 3)  # need to say 3 xi components expected
            element_identifier = element.getIdentifier()
            _, vagus_coordinates = marker_vagus_coordinates.evaluateReal(fieldcache, 3)
            if marker_host_coordinates:
                _, host_coordinates = marker_host_coordinates.evaluateReal(fieldcache, 3)
                data = (name, (element_identifier, xi), vagus_coordinates, host_coordinates)
            else:
                data = (name, (element_identifier, xi), vagus_coordinates)
            marker_data.append(data)
            node = nodeiterator.next()
        # by deleting the temporary field here, no clients will be notified of any change
        del marker_host_coordinates
    return marker_data
