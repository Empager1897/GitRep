import numpy as np
import math
import scipy

class Solution:
    def __init__(self, x0 = (1.0, 1.0, 1.0), default = True, eps = 1e-10):
        self._x0 = list(x0)
        self.__eps = eps if default else float(input("Enter eps (precision) : "))
        self._matrixA = np.array([
        [4.0, 1.0, 1.0],
        [1.0, 3.0, 0.0],
        [1.0, 0.0, 2.0],
    ]) if default else self.__generate_matrixA("Generating matrix A")
        self._matrixB = self.__generate_matrixB() if self._matrixA is not None else print("Matrix A is not suitable")

    def power_method(self)->tuple:
        lam_max_B = self.__lam_max()
        lam_min_A = self.__lam_min(self._matrixA, lam_max_B)
        lam_max_A = self.__lam_max(self._matrixA)
        return (lam_max_A, lam_min_A)

    def __lam_max(self, matrixB = None):
        lam = None
        lam_next = None
        lam_max_B = None
        counter = 0
        matrixB = matrixB if matrixB is not None else self._matrixB
        x = self._x0
        x_prev = self._x0
        while True :
            if counter > 1 and abs(lam_next - lam) < self.__eps :
                lam_max_B = lam_next
                print(f"Solution found after {counter} iterations")
                break
            x_prev = x
            x = matrixB @ x_prev
            lam = lam_next if counter > 0 else None
            lam_next = x[ 0 ] / x_prev[ 0 ]
            x_norm = self.__norm_vector(x)
            x = x / x_norm
            counter += 1
            print(f"Iteration {counter}: x1 = {x[ 0 ]}, x2 = {x[ 1 ]}, x3 = {x[ 2 ]} \n")
        return lam_max_B

    def __lam_min(self, matrix, lam_max):
        lam_min = self._matrix_infnorm(matrix) - lam_max
        return lam_min

    def __generate_matrixA(self, message):
        print(message)
        rows = int(input("Enter number of rows: "))
        cols = int(input("Enter number of columns: "))
        matrix = np.empty((rows,cols), dtype = int)
        check = True
        while check:
            for i in range(rows):
                for j in range(cols):
                    matrix[i][j] = float(input(f"Enter number a{i}{j} "))

            if np.array_equal(matrix, matrix.T):
                print("condition 1 of convergence satisfied (A_T = A)")
                check = False
            elif matrix.shape == (3, 1):
                check = False

        print("\n")
        return matrix

    def __generate_matrixB(self):
        print("Generating matrix B ")
        B = self._matrix_infnorm(self._matrixA) * np.eye(3) - self._matrixA
        print(f"B = {B} \n")
        return B

    def __norm_vector(self, vector):
        norm = 0.0
        for i in range(len(vector)):
            norm += vector[i] ** 2
        norm = math.sqrt(norm)
        return norm

    def _matrix_infnorm(self, matrix)->float:
        sums = []
        n = matrix.ndim
        for i in range(n):
            sum = 0.0
            for j in range(n):
                sum += abs(matrix[i,j])
            sums.append(sum)
        infnorm = max(sums)
        return infnorm

    def __private_method(self):
        pass


x0 = [1.0, 1.0, 1.0]
eps = 1e-10
ex1 = Solution(x0, True, eps)
lam_max_A, lam_min_A = ex1.power_method()
print(f"lam_max_A = {lam_max_A}, lam_min_A = {lam_min_A}")