import numpy as np

def format_complex(z):
    if z.imag == 0:
        return str(z.real)
    else:
        return str(z)

def matrix(n):
    filename=str(n)+'.csv'
    DirA="A/"
    DirB="B/"
    try:
        A = np.loadtxt(DirA+filename, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{DirA+filename}' が見つかりません。")
        return
    
    try:
        B = np.loadtxt(DirB+filename, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{DirB+filename}' が見つかりません。")
        return

    if(n==1):
        A=[[A]]

    det_A = np.linalg.det(A)
    B_col=B.reshape(-1)
    x = np.linalg.solve(A, B_col)
    eigvals,eigvecs=np.linalg.eig(A)
    
    with open("det/"+str(n)+'.csv', 'w', encoding='utf-8') as f:
        f.write(str(det_A))

    with open("x/"+str(n)+'.csv', 'w', encoding='utf-8') as f:
        f.write(','.join(map(str, x)) + '\n')

    with open("eigen/"+str(n)+'.csv', 'w', encoding='utf-8') as f:
        f.write("Index,Eigenvalue,")
        for i in range(n):
            f.write("Eigenvector_"+str(i))
            if i!=n-1:
                f.write(",")
        f.write("\n")
        for i in range(n):
            f.write(str(i)+","+format_complex(eigvals[i])+",")
            for j in range(n):
                f.write(format_complex(eigvecs[j][i]))
                if j!=n-1:
                    f.write(",")
            f.write("\n")
        


for i in range(100):
    matrix(i+1)