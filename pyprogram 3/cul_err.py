import numpy as np
import csv

#複素数の虚数部が0のとき、実数値の文字列で出力する
def format_complex(z):
    if z.imag == 0:
        return str(z.real)
    else:
        return str(z)

#eigenファイルを読み込み、eigvalとeigvecを返す
def read_eigen_csv(filename):
    eigenvalues = []
    eigenvectors = []

    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader) 

        for row in reader:
            # row[1] = 固有値, row[2:] = 固有ベクトル（n個）
            eigenvalues.append(complex(row[1]))  # 複素数対応
            vector = [complex(x) for x in row[2:]]
            eigenvectors.append(vector)

    eigenvalues = np.array(eigenvalues)
    eigenvectors = np.array(eigenvectors).T  # 列ベクトルとして扱うなら転置

    return eigenvalues, eigenvectors

#測定値と実際の値のノルム誤差を求める
def norm(eigvec,culeigvec):
    v1 = eigvec
    v2 = culeigvec

    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    if np.vdot(v1_norm, v2_norm).real < 0:
        v2_norm = -v2_norm
    error = np.linalg.norm(v1_norm - v2_norm)

    return error

#誤差を計算して出力
def cul_err(n):
    #csvファイルからデータを読み込む
    det_file="det/"+str(n)+".csv"
    try:
        det = np.loadtxt(det_file, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{det_file}' が見つかりません。")
        return
    
    culdet_file="cul_det/"+str(n)+".csv"
    try:
        culdet = np.loadtxt(culdet_file, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{culdet_file}' が見つかりません。")
        return

    x_file="x/"+str(n)+".csv"
    try:
        x = np.loadtxt(x_file, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{x_file}' が見つかりません。")
        return

    culx_file="cul_x/"+str(n)+".csv"
    try:
        culx = np.loadtxt(culx_file, delimiter=',', dtype=float)
    except FileNotFoundError:
        print(f"エラー: ファイル '{culx_file}' が見つかりません。")
        return

    culeigen_file="cul_eigen/"+str(n)+".csv"
    try:
        culeigvals,culeigvecs=read_eigen_csv(culeigen_file)
    except FileNotFoundError:
        print(f"エラー: ファイル '{culeigen_file}' が見つかりません。")
        return
    
    eigen_file="eigen/"+str(n)+".csv"
    try:
        eigvals,eigvecs=read_eigen_csv(eigen_file)
    except FileNotFoundError:
        print(f"エラー: ファイル '{eigen_file}' が見つかりません。")
        return

    if n==1:
        x=[x]
        culx=[culx]
        culeigvals=[culeigvals]
        culeigvecs=np.array(culeigvecs, dtype=complex).T
        eigvals=[eigvals]
        eigvecs=np.array(eigvecs, dtype=complex).T

    #誤差を求める
    EM_det_str=str(culdet/det)
    
    x.sort()
    culx.sort()
    EM_x_str=""
    EM_x=[]
    for i in range(n):
        EM_x.append(abs(1-culx[i]/x[i]))
        EM_x_str=EM_x_str+str(culx[i]/x[i])
        if i != n-1:
            EM_x_str=EM_x_str+","
    mean_EM_x=np.mean(EM_x)
    max_EM_x=np.max(EM_x)
    median_EM_x=np.median(EM_x)

    if n!=1:
        eigval_indices = np.argsort(eigvals)
        culeigval_indices = np.argsort(culeigvals)
        eigvals = eigvals[eigval_indices]
        eigvecs = eigvecs[:, eigval_indices]
        culeigvals = culeigvals[culeigval_indices]
        culeigvecs = culeigvecs[:, culeigval_indices]
    
    EM_eigvals_str=""
    EM_eigvals=[]
    for i in range(n):
        EM_eigvals.append(abs(1-culeigvals[i]/eigvals[i]))
        EM_eigvals_str=EM_eigvals_str+format_complex(culeigvals[i]/eigvals[i])
        if i != n-1:
            EM_eigvals_str=EM_eigvals_str+","
    
    mean_EM_eigvals=np.mean(EM_eigvals)
    max_EM_eigvals=np.max(EM_eigvals)
    median_EM_eigvals=np.median(EM_eigvals)
    
    EM_eigvecs_str=""
    EM_eigvecs=[]
    max_error=0
    for i in range(n):
        minimum_error=2
        for k in range(n):
            index=k
            if 0<=index<n:
                error=norm(eigvecs[:, index],culeigvecs[:, i])
                if(error<minimum_error):
                    minimum_error=error
        EM_eigvecs.append(minimum_error)
        EM_eigvecs_str=EM_eigvecs_str+str(minimum_error)
        if i != n-1:
            EM_eigvecs_str=EM_eigvecs_str+","
        if(max_error<minimum_error):
            max_error=minimum_error

    mean_EM_eigvecs=np.mean(EM_eigvecs)
    max_EM_eigvecs=np.max(EM_eigvecs)
    median_EM_eigvecs=np.median(EM_eigvecs)

    
    #結果を出力
    with open("err_measurement/"+str(n)+'.csv', 'w', encoding='utf-8') as f:
        f.write("det,"+EM_det_str+"\n")
        f.write("x:"+EM_x_str+"\n")
        f.write("max:"+str(max_EM_x)+", mean:"+str(mean_EM_x)+", median:"+str(median_EM_x)+"\n")
        f.write("eigvals,"+EM_eigvals_str+"\n")
        f.write("max:"+str(max_EM_eigvals)+", mean:"+str(mean_EM_eigvals)+", median:"+str(median_EM_eigvals)+"\n")
        f.write("eigvecs,"+EM_eigvecs_str+"\n")
        f.write("max:"+str(max_EM_eigvecs)+", mean:"+str(mean_EM_eigvecs)+", median:"+str(median_EM_eigvecs)+"\n")

    
for i in range(100):
    cul_err(i+1)