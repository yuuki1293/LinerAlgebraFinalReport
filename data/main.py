import numpy as np

def matrix(n):
    filename=str(n)+'.csv'
    DirA="A/"
    DirB="B/"
    DirX="x/"
    try:
        A = np.loadtxt(DirA+filename, delimiter=',', dtype=int)
    except FileNotFoundError:
        print(f"エラー: ファイル '{DirA+filename}' が見つかりません。")
        return
    
    try:
        B = np.loadtxt(DirB+filename, delimiter=',', dtype=int)
    except FileNotFoundError:
        print(f"エラー: ファイル '{DirB+filename}' が見つかりません。")
        return

    #try:
    #    x = np.loadtxt(DirX+filename, delimiter=',', dtype=int)
    #except FileNotFoundError:
    #    print(f"エラー: ファイル '{DirX+filename}' が見つかりません。")
    #    return

    det_A = np.linalg.det(A)
    B_col=B.reshape(-1)
    x = np.linalg.solve(A, B_col)
    eigvals,eigvecs=np.linalg.eig(A)

    for i in range(n):
        for j in range(n):
            print(f"{A[i][j]:>3}, ",end="")
        print()
    print()

    print("行列式:",det_A)
    print("連立方程式の解:",x)
    print("固有値:",eigvals)
    print("固有ベクトル:",eigvecs)


matrix(4)