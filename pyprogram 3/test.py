import os

# 拡張子の変換：.txt → .md
dir="/mnt/c/VScode/応用線形代数レポート/pyprogram/cul_det"
for filename in os.listdir(dir):
    new_name = f"{filename}.csv"
    new_full_path = os.path.join(dir, new_name)
    full_path = os.path.join(dir, filename)
    os.rename(full_path, new_full_path)