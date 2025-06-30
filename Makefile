# 応用線形代数 最終レポート Makefile

# コンパイラ設定
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -I./include
DEBUGFLAGS = -std=c++20 -Wall -Wextra -g -O0 -I./include -DDEBUG_MODE
LDFLAGS =

# ディレクトリ設定
SRCDIR = src
BUILDDIR = build
INCLUDEDIR = include

# ソースファイル
SOURCES = $(SRCDIR)/linear_algebra.cpp $(SRCDIR)/main.cpp
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# メインプログラム
MAIN_TARGET = $(BUILDDIR)/main

# デフォルトターゲット
all: $(MAIN_TARGET)

# デバッグビルド
debug: CXXFLAGS = $(DEBUGFLAGS)
debug: $(MAIN_TARGET)

# デバッグ実行（強制的にデバッグビルド）
debug-run: CXXFLAGS = $(DEBUGFLAGS)
debug-run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

# オブジェクトファイルの作成
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# メインプログラムの作成
$(MAIN_TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# メインプログラムの実行
run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

# グラフ生成
plot:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu
	gnuplot gnu/plot_eigenvalue_times.gnu
	gnuplot gnu/plot_linear_solver_times.gnu

# 個別計算時間のグラフ生成
plot-det:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu

plot-eigen:
	@mkdir -p figures
	gnuplot gnu/plot_eigenvalue_times.gnu

plot-linear:
	@mkdir -p figures
	gnuplot gnu/plot_linear_solver_times.gnu

# 全グラフ生成
plot-all:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu
	gnuplot gnu/plot_eigenvalue_times.gnu
	gnuplot gnu/plot_linear_solver_times.gnu

# クリーン
clean:
	rm -rf $(BUILDDIR)/*
	rm -rf data/*

# ヘルプ
help:
	@echo "利用可能なターゲット:"
	@echo "  all       - メインプログラムをビルド"
	@echo "  debug     - デバッグモードでビルド（n=10最大）"
	@echo "  run       - メインプログラムを実行"
	@echo "  debug-run - デバッグモードで実行"
	@echo "  plot      - 計算時間グラフを生成"
	@echo "  plot-det  - 行列式計算時間グラフを生成"
	@echo "  plot-eigen- 固有値・固有ベクトル計算時間グラフを生成"
	@echo "  plot-linear- 線形方程式解法時間グラフを生成"
	@echo "  plot-all  - 全グラフを生成"
	@echo "  clean     - ビルドファイルを削除"
	@echo "  help      - このヘルプを表示"

.PHONY: all debug run debug-run clean help
