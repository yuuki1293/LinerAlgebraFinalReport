# 固有値・固有ベクトル計算時間のプロット
set terminal pngcairo enhanced font "Arial,12" size 600,400
set output "../figures/eigen_times.png"
set datafile separator ","

# ラベル
set xlabel "行列サイズ" font "Arial,14"
set ylabel "計算時間 (ms)" font "Arial,14"

# グリッド
set grid

# データのプロット（ヘッダー行をスキップ）
plot "../data/eigen_times.csv" using 1:2 skip 1 with points pointtype 7 pointsize 0.7 title ""

# 対数スケール版
set output "../figures/eigen_times_log.png"
set logscale y
plot "../data/eigen_times.csv" using 1:2 skip 1 with points pointtype 7 pointsize 0.7 title ""
