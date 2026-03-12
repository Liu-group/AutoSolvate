## @file __main__.py
#  Gateway script to rest of program
#
# !/usr/bin/env python

import json
import os
import sys


def _detect_env_prefix() -> str | None:
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        return conda_prefix
    exe = sys.executable
    if exe:
        prefix = os.path.dirname(os.path.dirname(exe))
        if os.path.isdir(prefix):
            return prefix
    return None


def _build_mcp_server_config(server_name: str) -> dict:
    prefix = _detect_env_prefix()
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome and prefix:
        amber_bin = os.path.join(prefix, "bin")
        if os.path.exists(os.path.join(amber_bin, "antechamber")) or os.path.exists(os.path.join(amber_bin, "tleap")):
            amberhome = prefix
    env = {
        "PYTHONUNBUFFERED": "1",
    }
    if amberhome:
        env["AMBERHOME"] = amberhome
    if prefix:
        env["PATH"] = f"{prefix}/bin:${{env:PATH}}"
        env["LD_LIBRARY_PATH"] = f"{prefix}/lib:${{env:LD_LIBRARY_PATH}}"

    warnings = []
    if not prefix:
        warnings.append("Unable to detect environment prefix; PATH/LD_LIBRARY_PATH may be incomplete.")
    if not amberhome:
        warnings.append("AMBERHOME not detected; AmberTools-dependent workflows may fail.")
    if prefix:
        bin_dir = os.path.join(prefix, "bin")
        expected_bins = ["antechamber", "parmchk2", "tleap", "packmol", "obabel"]
        missing = [b for b in expected_bins if not os.path.exists(os.path.join(bin_dir, b))]
        if missing:
            warnings.append(f"Missing executables in {bin_dir}: {', '.join(missing)}")
    for message in warnings:
        print(f"Warning: {message}", file=sys.stderr)

    return {
        server_name: {
            "type": "stdio",
            "command": sys.executable,
            "args": ["-m", "autosolvate", "mcp_server", "--transport", "stdio"],
            "env": env,
        }
    }


def _emit_mcp_json(args: list[str]) -> None:
    server_name = "autosolvate-local"
    emit_full = False
    indent = 2
    i = 0
    while i < len(args):
        if args[i] == "--server-name" and i + 1 < len(args):
            server_name = args[i + 1]
            i += 2
            continue
        if args[i] == "--full":
            emit_full = True
            i += 1
            continue
        if args[i] == "--indent" and i + 1 < len(args):
            indent = int(args[i + 1])
            i += 2
            continue
        if args[i] in {"-h", "--help"}:
            print("Usage: autosolvate mcp_json [--full] [--server-name NAME] [--indent N]")
            print("  --full              emit a full mcp.json with a servers block")
            print("  --server-name NAME  MCP server name (default: autosolvate-local)")
            print("  --indent N          JSON indentation (default: 2)")
            return
        i += 1

    server_block = _build_mcp_server_config(server_name)
    payload = {"servers": server_block, "inputs": []} if emit_full else server_block[server_name]
    print(json.dumps(payload, indent=indent))

## Main function
#  @param args Argument namespace
def main(args=None):
    if args is None:
        args = sys.argv[1:]
    ### run with gui ###
    if len(args) == 0:
        print('AutoSolvate is starting with graphical user interface!')
        ### create main application
        from autosolvate.GUI.tk_autosolvate import autosolvateGUI, cleanUp, Tk
        window = Tk()
        my_gui = autosolvateGUI(window)
        window.mainloop()
        cleanUp()
    elif args[0] == '-h' or args[0] == '--help':
        print('Usage: autosolvate [OPTIONS]')
        print('  [NO OPTION]                launches GUI')
        print('  boxgen [OPTIONS]           generate initial structure')
        print('  boxgen_metal[OPTIONS]      generate initial structure and forcefield for organometallic compounds')
        print('  boxgen_multicomponent [OPTIONS] generate initial structure for multicomponent systems')
        print('  mdrun [OPTIONS]            automated QM/MM trajectory generation')
        print('  clustergen [OPTIONS]       extract microsolvated clusters ')
        print('  mcp_server [OPTIONS]       run FastMCP server for agent integration')
        print('  mcp_json [OPTIONS]         emit MCP server JSON for editor configs')
        print('  -h, --help                 short usage description')
        print()
        print('All options described at autosolvate.readthedocs.io')
        exit()
 
    elif args[0] == 'boxgen':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters.')
        from autosolvate.autosolvate import startboxgen
        startboxgen(args[1:])
    elif args[0] == 'boxgen_metal':
        from autosolvate.FFmetalcomplex import startFFgen
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters for organometallic compounds.')
        startFFgen(args[1:])
    elif args[0] == 'mdrun':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to automatically run MD simulations of solvated structure.')
        from autosolvate.generatetrajs import startmd
        startmd(args[1:])
    elif args[0] == 'clustergen':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to extract microsolvated clusters from MD trajectories with solvent box.')
        from autosolvate.clustergen import startclustergen
        startclustergen(args[1:])
    elif args[0] == 'boxgen_multicomponent':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters for multicomponent systems.')
        from autosolvate.multicomponent import startmulticomponent
        startmulticomponent(args[1:]) 
    elif args[0] == 'boxgen_interactive':
        print('AutoSolvate interactive input generator!')
        from autosolvate.cli.boxgen_interactive import main as interactive_main
        interactive_main(args[1:])
    elif args[0] == 'mcp_server':
        from autosolvate.fastmcp_server import main as mcp_main
        mcp_main(args[1:])
    elif args[0] == 'mcp_json':
        _emit_mcp_json(args[1:])
    else:
        print('Invalid syntax for AutoSolvate command line interface.')
        print('Please run \'autosolvate -h\' to check out the basic usage.')
        exit()

if __name__ == '__main__':
    main()
