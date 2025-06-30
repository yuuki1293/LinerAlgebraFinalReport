# 応用線形代数 最終レポート Makefile

# コンパイラ設定
CXX = g++
CXXFLAGS = -std=c++26 -Wall -Wextra -O2 -I./include
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

# クリーン
clean:
	rm -rf $(BUILDDIR)/*

# ヘルプ
help:
	@echo "利用可能なターゲット:"
	@echo "  all     - メインプログラムをビルド"
	@echo "  run     - メインプログラムを実行"
	@echo "  clean   - ビルドファイルを削除"
	@echo "  help    - このヘルプを表示"

.PHONY: all run clean help
