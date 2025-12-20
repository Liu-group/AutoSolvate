from __future__ import annotations

import argparse
import json
import os

from autosolvate.input_wizard import run_wizard
from autosolvate.utils.interactive_input_utils import InputAbort


def main(argv=None):
    parser = argparse.ArgumentParser(description="Interactive wizard for AutoSolvate multicomponent inputs")
    parser.add_argument("--agent-mode", action="store_true", help="Disable free-form edit mode for agents")
    parser.add_argument("-o", "--output", type=str, default="generated_input.json", help="Output JSON path")
    args = parser.parse_args(argv)

    try:
        config = run_wizard(agent_mode=args.agent_mode)
    except InputAbort:
        print("Wizard aborted by user.")
        return 0
    out_path = args.output
    if not os.path.isabs(out_path):
        out_path = os.path.join(config.get("folder", os.getcwd()), out_path)
    with open(out_path, "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"Saved input to {out_path}")


if __name__ == "__main__":
    main()
