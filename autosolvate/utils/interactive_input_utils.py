"""Utilities for the interactive input wizard."""
from __future__ import annotations

import ast
import os
from typing import Any, Callable, Dict, Iterable, List, Sequence, Tuple

MAX_RETRIES = 5


class InputAbort(Exception):
    """Raised when user chooses to exit the wizard."""


def _validate_control_words(raw: str) -> str:
    val = raw.strip()
    if val.lower() == "exit":
        raise InputAbort()
    return val


def ask_value(prompt: str, parser: Callable[[str], Any], validator: Callable[[Any], bool], default: Any = None, allow_skip: bool = False) -> Any:
    for _ in range(MAX_RETRIES):
        raw = input(prompt)
        raw = _validate_control_words(raw)
        if raw.strip() == "" and default is not None:
            return default
        if allow_skip and raw.lower() == "skip":
            return default
        try:
            value = parser(raw)
            if validator(value):
                return value
        except InputAbort:
            raise
        except Exception:
            pass
        print("Invalid input, please try again.")
    raise InputAbort()


def ask_yes_no(prompt: str, default: bool | None = None) -> bool:
    def parse(raw: str) -> str:
        return _validate_control_words(raw).strip().lower()

    def valid(val: str) -> bool:
        return val in {"yes", "no", ""}

    for _ in range(MAX_RETRIES):
        raw = input(prompt)
        raw = _validate_control_words(raw)
        val = raw.strip().lower()
        if val == "" and default is not None:
            return default
        if val in {"yes", "no"}:
            return val == "yes"
        print("Please answer yes or no.")
    raise InputAbort()


def parse_cube_size(raw: str):
    parts = [p.strip() for p in raw.replace(";", ",").split(",") if p.strip()]
    if len(parts) == 1:
        val = float(parts[0])
        if val <= 0:
            raise ValueError
        return val
    if len(parts) == 3:
        vals = [float(p) for p in parts]
        if any(v <= 0 for v in vals):
            raise ValueError
        return vals
    raise ValueError


def parse_int(raw: str) -> int:
    return int(raw.strip())


def parse_float(raw: str) -> float:
    return float(raw.strip())


def basename_no_ext(path: str) -> str:
    base = os.path.basename(path)
    return os.path.splitext(base)[0]


def apply_config_edit(config: dict, command: str) -> bool:
    """Apply a minimal config[...]=... edit; returns True if applied."""
    tree = ast.parse(command, mode="exec")
    if len(tree.body) != 1 or not isinstance(tree.body[0], ast.Assign):
        return False
    assign = tree.body[0]
    if len(assign.targets) != 1:
        return False
    target = assign.targets[0]
    if not isinstance(target, ast.Subscript) or not isinstance(target.value, ast.Name) or target.value.id != "config":
        return False
    # collect key chain
    keys: List[Any] = []
    while isinstance(target, ast.Subscript):
        if isinstance(target.slice, ast.Index):  # type: ignore[attr-defined]
            key_node = target.slice.value
        else:
            key_node = target.slice
        if isinstance(key_node, (ast.Constant, ast.Str)):
            keys.append(key_node.value if hasattr(key_node, "value") else key_node.s)
        else:
            return False
        if isinstance(target.value, ast.Subscript):
            target = target.value
        else:
            break
    value = ast.literal_eval(assign.value)
    _set_nested(config, list(reversed(keys)), value)
    return True


def _set_nested(cfg: dict, keys: List[Any], value: Any):
    cur = cfg
    for k in keys[:-1]:
        if isinstance(k, int):
            while len(cur) <= k:
                cur.append({})
            if cur[k] is None:
                cur[k] = {}
            cur = cur[k]
        else:
            if k not in cur or cur[k] is None:
                cur[k] = {}
            cur = cur[k]
    last = keys[-1]
    if isinstance(last, int):
        while len(cur) <= last:
            cur.append(None)
        cur[last] = value
    else:
        cur[last] = value
