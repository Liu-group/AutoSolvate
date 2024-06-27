NameError: name '_read_utf8_with_fallback' is not defined



(test2) [pli@pascal run1]$ autosolvate boxgen_multicomponent
Traceback (most recent call last):
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 589, in _build_master
    ws.require(__requires__)
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 926, in require
    needed = self.resolve(parse_requirements(requirements))
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 787, in resolve
    dist = self._resolve_dist(
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 833, in _resolve_dist
    raise VersionConflict(dist, req).with_context(dependent_req)
pkg_resources.VersionConflict: (autosolvate 0.1.4+136.gda10105.dirty (/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/autosolvate-0.1.4+136.gda10105.dirty-py3.8.egg), Requirement.parse('autosolvate==0.1.4+154.gd24f652.dirty'))

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/pli/anaconda3/envs/test2/bin/autosolvate", line 33, in <module>
    sys.exit(load_entry_point('autosolvate==0.1.4+154.gd24f652.dirty', 'console_scripts', 'autosolvate')())
  File "/home/pli/anaconda3/envs/test2/bin/autosolvate", line 25, in importlib_load_entry_point
    return next(matches).load()
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/importlib/metadata.py", line 77, in load
    module = import_module(match.group('module'))
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/importlib/__init__.py", line 127, in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
  File "<frozen importlib._bootstrap>", line 1014, in _gcd_import
  File "<frozen importlib._bootstrap>", line 991, in _find_and_load
  File "<frozen importlib._bootstrap>", line 961, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 219, in _call_with_frames_removed
  File "<frozen importlib._bootstrap>", line 1014, in _gcd_import
  File "<frozen importlib._bootstrap>", line 991, in _find_and_load
  File "<frozen importlib._bootstrap>", line 975, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 671, in _load_unlocked
  File "<frozen importlib._bootstrap_external>", line 843, in exec_module
  File "<frozen importlib._bootstrap>", line 219, in _call_with_frames_removed
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/autosolvate-0.1.4+136.gda10105.dirty-py3.8.egg/autosolvate/__init__.py", line 7, in <module>
    from .autosolvate import *
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/autosolvate-0.1.4+136.gda10105.dirty-py3.8.egg/autosolvate/autosolvate.py", line 1, in <module>
    from .molecule import *
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/autosolvate-0.1.4+136.gda10105.dirty-py3.8.egg/autosolvate/molecule/__init__.py", line 1, in <module>
    from .molecule import *
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/autosolvate-0.1.4+136.gda10105.dirty-py3.8.egg/autosolvate/molecule/molecule.py", line 6, in <module>
    import pkg_resources
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 3283, in <module>
    def _initialize_master_working_set():
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 3266, in _call_aside
    f(*args, **kwargs)
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 3295, in _initialize_master_working_set
    working_set = _declare_state('object', 'working_set', WorkingSet._build_master())
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 591, in _build_master
    return cls._build_from_requirements(__requires__)
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 604, in _build_from_requirements
    dists = ws.resolve(reqs, Environment())
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 1014, in __init__
    self.scan(search_path)
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 1046, in scan
    for dist in find_distributions(item):
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 2091, in find_on_path
    yield from factory(fullpath)
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 2183, in resolve_egg_link
    return next(dist_groups, ())
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 2179, in <genexpr>
    resolved_paths = (
  File "/home/pli/anaconda3/envs/test2/lib/python3.8/site-packages/pkg_resources/__init__.py", line 2167, in non_empty_lines
    for line in _read_utf8_with_fallback(path).splitlines():
NameError: name '_read_utf8_with_fallback' is not defined








autosolvate mdrun will conflict with results generated by autosovlate boxgen_multicomponent







cant use nohup to run the following command 

autosolvate mdrun -f MYBOX -q 0 -u 1 -t 300 -p 1 -m 10000 -n 10000 -o 100 -s 100 -l 250 -r 



 