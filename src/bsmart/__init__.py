
import os
import sys



# --- METADATA ---

import importlib.metadata
__version__ = "PRIVATEVERSION"
try:
    # This must match the 'name' field in your pyproject.toml
    __version__ = importlib.metadata.version("bsmart")
except importlib.metadata.PackageNotFoundError:
    # This handles the case where the package is imported locally
    # but not actually installed (e.g., during some CI steps or local scripts)
    __version__ = "PRIVATEVERSION"



__build_date__ = "RELEASEDATE"




def get_bsmart_config_dir(create_dir=False):
    """
    Returns the path to the BSMArt configuration directory.
    Order of preference:
    1. Virtual Environment: sys.prefix/etc/BSMArt/
    2. User Config: XDG_CONFIG_HOME/BSMArt/ (e.g. ~/.config/BSMArt/)
    """
    
    # Check for virtual environment
    # sys.prefix != sys.base_prefix is the standard way to check for venv/conda in Python 3.3+
    in_venv = (sys.prefix != getattr(sys, "base_prefix", sys.prefix))
    
    if in_venv:
        # Use venv-local config
        bsmart_config_dir = os.path.join(sys.prefix, 'etc', 'BSMArt')
    else:
        # Use global user config
        config_home = os.environ.get('XDG_CONFIG_HOME')
        if not config_home:
            config_home = os.path.join(os.path.expanduser('~'), '.config')
        bsmart_config_dir = os.path.join(config_home, 'BSMArt')
    
    if create_dir and not os.path.exists(bsmart_config_dir):
        try:
            os.makedirs(bsmart_config_dir)
        except OSError:
            # Fallback or error logging could go here, but for now we let it fail or user handle it
            pass
            
    return bsmart_config_dir


def get_heptools_path_file(create_dir=False):
    config_dir=get_bsmart_config_dir(create_dir=create_dir)
    return os.path.join(config_dir, 'HEPtoolpaths.json')

# def get_heptools_path_file(create_dir=False):
#     """
#     Returns the path to the HEPtoolpaths.json file.
#     Order of preference:
#     1. Virtual Environment: sys.prefix/etc/BSMArt/HEPtoolpaths.json
#     2. User Config: XDG_CONFIG_HOME/BSMArt/HEPtoolpaths.json (e.g. ~/.config/BSMArt/HEPtoolpaths.json)
#     """
    
#     # Check for virtual environment
#     # sys.prefix != sys.base_prefix is the standard way to check for venv/conda in Python 3.3+
#     in_venv = (sys.prefix != getattr(sys, "base_prefix", sys.prefix))
    
#     if in_venv:
#         # Use venv-local config
#         bsmart_config_dir = os.path.join(sys.prefix, 'etc', 'BSMArt')
#     else:
#         # Use global user config
#         config_home = os.environ.get('XDG_CONFIG_HOME')
#         if not config_home:
#             config_home = os.path.join(os.path.expanduser('~'), '.config')
#         bsmart_config_dir = os.path.join(config_home, 'BSMArt')
    
#     if create_dir and not os.path.exists(bsmart_config_dir):
#         try:
#             os.makedirs(bsmart_config_dir)
#         except OSError:
#             # Fallback or error logging could go here, but for now we let it fail or user handle it
#             pass
            
#     return os.path.join(bsmart_config_dir, 'HEPtoolpaths.json')

