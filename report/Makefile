# TeXレポート用Makefile

# コンパイラ設定
LATEX = platex
DVIPDFMX = dvipdfmx
BIBTEX = pbibtex

# ファイル名
TARGET = main
SOURCE = $(TARGET).tex

# デフォルトターゲット
all: $(TARGET).pdf

# PDFファイルの生成
$(TARGET).pdf: $(TARGET).dvi
	$(DVIPDFMX) $(TARGET).dvi

# DVIファイルの生成
$(TARGET).dvi: $(SOURCE)
	$(LATEX) $(SOURCE)
	$(LATEX) $(SOURCE)  # 2回実行して相互参照を解決

# クリーンアップ
clean:
	rm -f *.aux *.log *.dvi *.pdf *.toc *.out *.bbl *.blg *.fdb_latexmk *.fls *.synctex.gz

# 完全クリーンアップ
distclean: clean
	rm -f $(TARGET).pdf

# ヘルプ
help:
	@echo "利用可能なコマンド:"
	@echo "  make        - PDFファイルを生成"
	@echo "  make clean  - 中間ファイルを削除"
	@echo "  make distclean - すべての生成ファイルを削除"
	@echo "  make help   - このヘルプを表示"

.PHONY: all clean distclean help
