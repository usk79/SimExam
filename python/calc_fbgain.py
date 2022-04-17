from control.matlab import *
import numpy as np

def check_ctrb(A, B):
    Uc = ctrb(A, B) # 可制御性行列の計算
    Nu = np.linalg.matrix_rank(Uc)  # Ucのランクを計算
    (N, N) = np.matrix(A).shape     # 正方行列Aのサイズ(N*N)
    # 可制御性の判別
    if Nu == N: return 0            # 可制御
    else: return -1                 # 可制御でない

def main():
    # RLC直列回路の状態空間表現
    r = 10.0
    l = 100e-3
    c = 100e-6

    A = np.array([[-r / l, -1 / (l * c)],
                [1.0   ,  0.0]])
    B = np.array([[1.0 / l], [0.0]])

    if check_ctrb(A, B) == -1 :
        print("uncontrolable") 
        exit

    print("A = ", A)
    print("B = ", B)

    poles = [-20 + 10j, -20 - 10j]

    F = place(A, B, poles)

    print(F)


if __name__ == "__main__":
  main()