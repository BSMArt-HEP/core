import pkgutil
import ast
import json
import importlib.util
import sys
from pathlib import Path
import os 
def _extract_metadata(source):
    """
    Parse the source code and return the Python value assigned to __meta__,
    using ast to avoid executing code. Returns None if not found.
    """
    try:
        # Parse module into an AST
        module = ast.parse(source)
    except SyntaxError as e:  
        return {"Syntax error" : str(e)}, None

    documentation =  ast.get_docstring(module)  # returns None if there is none

    # Look for top-level assignments to __meta__
    for node in module.body:
        # We only consider simple Assign nodes (not AugAssign, Not annotated assigns in 3.5)
        if isinstance(node, ast.Assign):
            # node.targets can be a list; check each
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__meta__":
                    # node.value is an AST expression for the RHS
                    try:
                        # Safely evaluate literal structures (dicts, lists, strings, numbers, booleans, None)
                        value = ast.literal_eval(node.value)
                        return value, documentation
                    except Exception:
                        # If literal_eval fails (e.g. non-literal code), return error marker
                        return {"error": "Non-literal __meta__ value; cannot safely evaluate."}, None
    return {}, documentation

def get_available_scans(search_target=None, with_meta:bool = False) -> dict:
    """
    Return a list of available scan modules and their metadata.
    search_target can be:
      - None: searches the current package (bsmart.scans)
      - A module object (e.g. bsmart.tools): searches that package
      - A dictionary path (str): searches that directory (e.g. ./Tools)
    """
    
    paths = None
    prefix = ""
    base_import = None
    
    # 1. Determine paths and prefix based on input type
    if search_target is None:
        # Default: default to this package (bsmart.scans)
        # We avoid using the global __path__ directly if possible, or use it if we are sure.
        # Ideally, we import ourselves to be sure?
        # But __path__ is available in scope.
        paths = __path__
        prefix = __name__ + "."
        base_import = None
            
    elif hasattr(search_target, '__path__'):
        # It's a package module (e.g. bsmart.tools)
        paths = search_target.__path__
        prefix = search_target.__name__ + "."
        base_import = search_target.__name__ # or None?

    elif isinstance(search_target, str):
        # Assume it's a directory path (legacy/local tools)
        if os.path.isdir(search_target):
            paths = [search_target]
            prefix = ""
            base_import = search_target
        else:
            # Maybe it's a package name string? Not handling for now unless requested.
            return {}
    else:
        # Unknown input type
        return {}

    scans = {}

    # 2. Walk packages
    for module_finder, name, ispkg in pkgutil.walk_packages(paths, prefix=prefix):
        if ispkg:
            continue

        # 3. Determine scan_id
        # internal scans (bsmart.scans.X) -> id = X
        # tools (bsmart.tools.Y) -> id = Y
        # local files (X) -> id = X
        
        if prefix and name.startswith(prefix):
            scan_id = name[len(prefix):]
        else:
            scan_id = name

        scan_entry = {"module": name, "import_path": base_import}
        
        # 4. Try to find path for display
        try:
             if hasattr(module_finder, 'path'):
                 # Construct valid path if possible, mostly for user info
                 scan_entry['path'] = str(Path(module_finder.path) / (name.split('.')[-1] + ".py"))
             else:
                 scan_entry['path'] = "Zip/Unknown"
        except:
             scan_entry['path'] = "Unknown"

        # 5. Extract metadata if requested
        if with_meta:
            source = None
            try:
                # Try to get source using the finder/loader
                if hasattr(module_finder, 'find_spec'):
                    spec = module_finder.find_spec(name)
                    if spec and spec.loader:
                        source = spec.loader.get_source(name)
                elif hasattr(module_finder, 'find_module'):
                    loader = module_finder.find_module(name)
                    if loader:
                        source = loader.get_source(name)
            except Exception:
                pass

            if source:
                meta, docs = _extract_metadata(source)
                scan_entry['meta'] = meta
                if docs is not None:
                    scan_entry['documentation'] = docs

        scans[scan_id] = scan_entry
    
    return scans
