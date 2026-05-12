#!/usr/bin/env python3

import argparse
import os
import sys
import collections
import importlib.util
import json
import ast

# Package imports
import bsmart
import bsmart.tools
from bsmart.scans import get_available_scans
from bsmart.scripts.bsmart_run import parse_input

def check_package(package_name):
    """
    Check if a package is installed.
    Returns (True, None) if installed.
    Returns (False, ErrorMessage) if not installed.
    """
    try:
        spec = importlib.util.find_spec(package_name)
        if spec is None:
             return False, f"Package '{package_name}' not found."
        return True, None
    except ImportError:
         return False, f"Package '{package_name}' could not be imported."
    except Exception as e:
         return False, f"Error checking package '{package_name}': {e}"

def parse_meta(meta_data):
    """
    Parse the metadata which can be a dict or a string.
    If it's a string, we assume it's the body of a JSON object (attributes not enclosed in {}) 
    or a python dict string.
    """
    if isinstance(meta_data, dict):
        return meta_data
    
    if isinstance(meta_data, str):
        # Try to parse as JSON first (wrapping in braces)
        try:
            # SModelS example: "name": "...", ...
            # We wrap it in {}
            json_str = "{" + meta_data + "}"
            return json.loads(json_str)
        except json.JSONDecodeError:
            pass
        
        # Try evaluating as python literal
        try:
            return ast.literal_eval("{" + meta_data + "}")
        except Exception:
            pass
            
        # Try evaluating as python literal without adding braces (maybe it already has them?)
        try:
            return ast.literal_eval(meta_data)
        except Exception:
            pass

    return {}

def main():
    parser = argparse.ArgumentParser(
        description='Check if required packages are installed for BSMArt scans and tools.')
    
    # Mode 1: Input file
    parser.add_argument('inputfile', metavar='File', type=str, nargs='?',
                        help='Input JSON file to parse for scan and tools', default=None)
    
    # Mode 2: Command line args
    parser.add_argument("--scan", help="Check packages for specific scan(s)",
                        nargs='+', action="append", default=[])
    parser.add_argument("--tool", help="Check packages for specific tool(s)",
                        nargs='+', action="append", default=[])
    parser.add_argument("--install", help="Install missing packages",
                        action="store_true")
    
    args = parser.parse_args()

    # Flatten lists of lists from action='append' + nargs='+'
    scans_to_check = set()
    for scan_list in args.scan:
        scans_to_check.update(scan_list)
        
    tools_to_check = set()
    for tool_list in args.tool:
        tools_to_check.update(tool_list)
    
    # If input file is provided, add scans/tools from it
    if args.inputfile:
        if not os.path.isfile(args.inputfile):
            print(f"Error: Input file '{args.inputfile}' not found.")
            sys.exit(1)
        
        try:
            inputs = parse_input(args.inputfile)
            
            # Scans
            if 'Setup' in inputs and 'Type' in inputs['Setup']:
                scans_to_check.add(inputs['Setup']['Type'])
            
            # Tools
            if 'Codes' in inputs:
                for tool in inputs['Codes']:
                    tools_to_check.add(tool)
                    
        except Exception as e:
            print(f"Error parsing input file: {e}")
            sys.exit(1)

    if not scans_to_check and not tools_to_check:
        print("No scans or tools specified. Usage:")
        print("  BSMArt-CheckPackages [--install] <inputfile.json>")
        print("  BSMArt-CheckPackages [--install] --scan <ScanName> [ScanName2 ...] --tool <ToolName> [ToolName2 ...]")
        sys.exit(0)

    # Gather available scans and tools
    # Note: get_available_scans returns both scans and tools depending on where we look
    
    # Local paths
    cwd = os.getcwd()
    local_scanpath = os.path.join(cwd, 'Scans')
    local_toolpath = os.path.join(cwd, 'Tools')
    
    # 1. Available Scans (package + local)
    available_scans = get_available_scans(None, with_meta=True)
    available_scans.update(get_available_scans(local_scanpath, with_meta=True))
    
    # 2. Available Tools (package + local)
    # Package tools
    main_tooldir = bsmart.tools
    available_tools = get_available_scans(main_tooldir, with_meta=True)
    available_tools.update(get_available_scans(local_toolpath, with_meta=True))

    missing_packages = set()
    all_ok = True
    
    print("Checking packages...")
    print("-" * 40)

    # Check Scans
    for scan_name in scans_to_check:
        if scan_name not in available_scans:
            print(f"Scan '{scan_name}': [NOT FOUND]")
            all_ok = False
            continue
            
        entry = available_scans[scan_name]
        meta = parse_meta(entry.get('meta', {}))
        
        requires = meta.get('requires', [])
        if not requires:
            print(f"Scan '{scan_name}': No requirements found.")
        else:
            print(f"Scan '{scan_name}' requirements:")
            for req in requires:
                ok, msg = check_package(req)
                status = "OK" if ok else f"MISSING ({msg})"
                print(f"  - {req}: [{status}]")
                if not ok:
                    all_ok = False
                    missing_packages.add(req)

    # Check Tools
    for tool_name in tools_to_check:
        if tool_name not in available_tools:
            print(f"Tool '{tool_name}': [NOT FOUND]")
            all_ok = False
            continue
            
        entry = available_tools[tool_name]
        meta = parse_meta(entry.get('meta', {}))
        
        requires = meta.get('requires', [])
        if not requires:
            print(f"Tool '{tool_name}': No requirements.")
        else:
            print(f"Tool '{tool_name}' requirements:")
            for req in requires:
                ok, msg = check_package(req)
                status = "OK" if ok else f"MISSING ({msg})"
                print(f"  - {req}: [{status}]")
                if not ok:
                    all_ok = False
                    missing_packages.add(req)
                    
    print("-" * 40)

    # Install missing packages if requested
    if args.install and missing_packages:
        print(f"Installing missing packages: {', '.join(missing_packages)}...")
        print("-" * 40)
        import subprocess
        
        for pkg in missing_packages:
            print(f"Installing '{pkg}'...")
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
                print(f"Successfully installed '{pkg}'.")
            except subprocess.CalledProcessError as e:
                print(f"Failed to install '{pkg}': {e}")
            except Exception as e:
                print(f"Error installing '{pkg}': {e}")
        
        print("-" * 40)
        print("Re-checking packages after installation...")
        print("-" * 40)
        
        # Re-check logic
        all_ok = True 
        # We re-check everything to be thorough/consistent, but simplest is just checking the missing ones
        # However, for output consistency, let's re-verify the originally missing ones.
        for pkg in missing_packages:
             ok, msg = check_package(pkg)
             status = "OK" if ok else f"MISSING ({msg})"
             print(f"  - {pkg}: [{status}]")
             if not ok:
                 all_ok = False

        print("-" * 40)

    if all_ok:
        print("All checks passed!")
        sys.exit(0)
    else:
        print("Some checks failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
