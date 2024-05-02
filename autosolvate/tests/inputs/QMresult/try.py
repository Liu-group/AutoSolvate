import pkgutil
import os

# 获取资源的二进制内容
data = pkgutil.get_data('autosolvate', os.path.join('data', solvPrefix, solvent_pdb))

# 如果你需要将这些数据写入临时文件以获取路径
import tempfile

with tempfile.NamedTemporaryFile('wb', delete=False) as tmp_file:
    tmp_file.write(data)
    solvent_pdb_origin = tmp_file.name
print(solvent_pdb_origin)
