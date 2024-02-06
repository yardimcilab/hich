import itertools

def fmtl(lst, join_char = " ", flag = ""):
    return f"{flag} {join_char.join(lst)}" if lst else ""

def _constrain_to(yaml, keys):
    constraints = []
    next_key = keys[0]
    keys_remain = keys[1:]
    assert next_key in yaml, f"In setting constraints, {next_key} not in {yaml}. Check config.yaml."
    if keys_remain:
        for key, val in yaml[next_key].items():
            constraints += _constrain_to(val, keys_remain)
        return constraints
    else:
        if isinstance(yaml[next_key], dict):
            return list(yaml[next_key].keys())
        elif isinstance(yaml[next_key], list):
            return yaml[next_key]
        else:
            return [yaml[next_key]]

def constrain_to(yaml, keys):
    assert keys, f"In setting constraints, no constraint-setting keys provided. Check call to constrain_to."
    assert yaml, f"In setting constraints, no constraint-setting keys provided. Check call to constrain_to."
    constraints = _constrain_to(yaml, keys)
    constraints_unique = list(set(constraints))
    return '|'.join(constraints_unique)

def root_to_leaf_paths(tree, path = []):
    """
    Generates all root-to-leaf paths in a given tree structure.

    :param tree: A dictionary representing the tree structure.
    :param path: A list representing the current path. This should not be set by the user.
    :return: A list of lists, where each inner list is a path from the root to a leaf in the tree.
    """
    result = []
    if not tree:
        return []

    for node, subtree in tree.items():
        current_path = path + [node]
        if isinstance(subtree, dict):
            result.extend(root_to_leaf_paths(subtree, current_path))
        else:
            result.append(current_path + [subtree])
    return result

def predecessor_successor(root_leaf_path, predecessors):
    """
    Creates a dictionary mapping each predecessor in the given path to its successor.

    :param root_leaf_path: A list representing a path from root to leaf in a tree.
    :param predecessors: A list of nodes for which successors are to be found.
    :return: A dictionary where each key is a predecessor and its value is the respective successor.
             Returns an empty dictionary if any predecessor is not found in the path.
    """
    mapping = {}
    for predecessor in predecessors:
        if predecessor not in root_leaf_path[:-1]:
            return {}
        mapping[predecessor] = successor(root_leaf_path, predecessor)
    return mapping

def successor(lst: list, predecessor):
    """
    Finds and returns the successor of a given element in a list.

    :param lst: The list in which to search for the successor.
    :param predecessor: The element whose successor is to be found.
    :return: The successor of the given element in the list.
    :raises ValueError: If the predecessor is not in the list or is the last element.
    """
    try:
        index = lst.index(predecessor)
        if index == len(lst) - 1:
            raise ValueError("The given predecessor is the last element, hence it has no successor.")
        return lst[index + 1]
    except ValueError:
        raise ValueError(f"Predecessor '{predecessor}' not found in the list.")
        
# Reimplementation of expand
def _expand(template, **kwargs):
    """
    Generates formatted strings by expanding the template with all combinations of keyword arguments.
    Reimplementation of Snakemake's 'expand' function to permit use and unit-testing outside a Snakemake context
    
    :param template: A string template with placeholders for formatting.
    :param kwargs: Keyword arguments where each key has a list of values to be combined.
    :return: A list of formatted strings.
    """
    keys, values = zip(*kwargs.items())
    combinations = itertools.product(*values)
    return [template.format(**dict(zip(keys, combo))) for combo in combinations]


def tree_expand(template, tree, tree_template_kw, **other_template_kwargs):
    """
    Expands a template string based on the paths in a tree structure and additional keyword arguments.
    
    :param template: A string template for formatting.
    :param tree: A tree structure, such as YAML.
    :param tree_template_kw: A list of keyword mappings to be used in the template.
                             Can be bare strings or tuples. If a tuple (TREE, TEMPLATE)
                             is used, TREE is the tree/YAML keyword and TEMPLATE is
                             the corresponding wildcard in the template string. If a bare
                             string is submitted, the same keyword will be used for both.
    :param other_template_kwargs: Other keyword arguments for the template.
    :return: A list of uniquely expanded strings.

    """
    if not tree_template_kw:
        return _expand(template, **other_template_kwargs)

    tree_template_kw = [(ttkw, ttkw) if isinstance(ttkw, str) else ttkw for ttkw in tree_template_kw]
    tree_kw, template_kw = zip(*tree_template_kw)
    results = set()

    for path in root_to_leaf_paths(tree):
        key_val = predecessor_successor(path, tree_kw)
        if key_val:
            new_kwargs = other_template_kwargs.copy()

            # Map template keywords to tree values
            for i, kw in enumerate(template_kw):
                new_kwargs[kw] = key_val[tree_kw[i]]

            # Convert kwarg bare strings to a single item list for expand to function properly
            new_kwargs = {k: [v] if not isinstance(v, list) else v for k, v in new_kwargs.items()}
            results.update(_expand(template, **new_kwargs))
    return list(results)