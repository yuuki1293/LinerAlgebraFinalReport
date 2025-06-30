# 固有値・固有ベクトル計算時間のプロット
set terminal pngcairo enhanced font "Arial,12" size 800,600
set output "figures/eigenvalue_times.png"
set datafile separator ","

# ラベル
set xlabel "行列サイズ" font "Arial,14"
set ylabel "計算時間 (ms)" font "Arial,14"

# グリッド
set grid

# データのプロット（ヘッダー行をスキップ）
plot "data/eigenvalue_times.csv" using 1:2 skip 1 with points pointtype 7 pointsize 0.7 title "固有値・固有ベクトル計算時間"

# 対数スケール版
set output "figures/eigenvalue_times_log.png"
set logscale y
plot "data/eigenvalue_times.csv" using 1:2 skip 1 with points pointtype 7 pointsize 0.7 title "固有値・固有ベクトル計算時間"
