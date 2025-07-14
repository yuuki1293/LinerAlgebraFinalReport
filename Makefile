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
TEST_SOURCES = $(SRCDIR)/linear_algebra.cpp $(SRCDIR)/test_functions.cpp
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)
DEBUG_OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%-debug.o)
TEST_OBJECTS = $(TEST_SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/test-%.o)

# メインプログラム
MAIN_TARGET = $(BUILDDIR)/main
DEBUG_TARGET = $(BUILDDIR)/main-debug
TEST_TARGET = $(BUILDDIR)/test

# デフォルトターゲット
all: $(MAIN_TARGET)

# デバッグビルド
debug: CXXFLAGS = $(DEBUGFLAGS)
debug: $(DEBUG_TARGET)

# デバッグ実行（強制的にデバッグビルド）
debug-run: CXXFLAGS = $(DEBUGFLAGS)
debug-run: $(DEBUG_TARGET)
	./$(DEBUG_TARGET)

# テストビルド
test: CXXFLAGS = $(DEBUGFLAGS)
test: $(TEST_TARGET)

# テスト実行
test-run: CXXFLAGS = $(DEBUGFLAGS)
test-run: $(TEST_TARGET)
	./$(TEST_TARGET)

# オブジェクトファイルの作成
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# デバッグオブジェクトファイルの作成
$(BUILDDIR)/%-debug.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# テストオブジェクトファイルの作成
$(BUILDDIR)/test-%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# メインプログラムの作成
$(MAIN_TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# デバッグメインプログラムの作成
$(DEBUG_TARGET): $(DEBUG_OBJECTS)
	$(CXX) $(DEBUG_OBJECTS) -o $@ $(LDFLAGS)

# テストプログラムの作成
$(TEST_TARGET): $(TEST_OBJECTS)
	$(CXX) $(TEST_OBJECTS) -o $@ $(LDFLAGS)

# メインプログラムの実行
run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

# 固有値・固有ベクトルグラフ生成
eigen_plots:
	@mkdir -p figures
	gnuplot gnu/plot_eigen_times.gnu

# 計算時間のグラフ生成
plot:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu
	gnuplot gnu/plot_linear_solver_times.gnu
	$(MAKE) eigen_plots

# 個別計算時間のグラフ生成
plot-det:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu

plot-eigen:
	@mkdir -p figures
	gnuplot gnu/plot_eigen_times.gnu

plot-linear:
	@mkdir -p figures
	gnuplot gnu/plot_linear_solver_times.gnu

# 全グラフ生成
plot-all:
	@mkdir -p figures
	gnuplot gnu/plot_determinant_times.gnu
	gnuplot gnu/plot_linear_solver_times.gnu
	$(MAKE) eigen_plots

# クリーン
clean:
	rm -rf $(BUILDDIR)/*
	rm -rf data/*

# ヘルプ
help:
	@echo "利用可能なターゲット:"
	@echo "  all       - メインプログラムをビルド"
	@echo "  debug     - デバッグモードでビルド（n=10最大）"
	@echo "  test      - 単体テストをビルド"
	@echo "  test-run  - 単体テストを実行"
	@echo "  run       - メインプログラムを実行"
	@echo "  debug-run - デバッグモードで実行"
	@echo "  plot      - 計算時間グラフを生成"
	@echo "  plot-det  - 行列式計算時間グラフを生成"
	@echo "  plot-eigen- 固有値・固有ベクトル計算時間グラフを生成"
	@echo "  plot-linear- 線形方程式解法時間グラフを生成"
	@echo "  plot-all  - 全グラフを生成"
	@echo "  clean     - ビルドファイルを削除"
	@echo "  help      - このヘルプを表示"

.PHONY: all debug test test-run run debug-run clean help
