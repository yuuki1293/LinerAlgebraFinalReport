# 条件数のプロット
set terminal pngcairo enhanced font "Arial,12" size 800,600
set output "../figures/condition_numbers.png"
set datafile separator ","

# ラベル
set xlabel "行列サイズ" font "Arial,14"
set ylabel "条件数" font "Arial,14"

# グリッド
set grid

# データのプロット（ヘッダー行をスキップ）
plot "../condition_numbers.csv" using 1:2 skip 1 with points pointtype 7 pointsize 0.7 title ""
