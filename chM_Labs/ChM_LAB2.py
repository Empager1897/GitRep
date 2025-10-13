import numpy as np
import pandas as pd
import sys
import math
import pprint

np.set_printoptions(threshold = sys.maxsize, precision = 2, suppress = True, linewidth = 100)

class Solution:
    def __init__(self):
        self.matrixA = self.generate_matrix("Printing A matrix ")
        self.rightVector = self.generate_matrix("Printing right vector b ")
        self.matrixS = np.zeros((3 , 3) , dtype = float)
        self.matrixD = np.zeros((3 , 3) , dtype = float)

    def generate_matrix(self, message):
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

    def sumforSD(self, i, j):
        sumdiag = 0
        sumnotdiag = 0
        for p in range(i):
            sumdiag += self.matrixS[p, i] ** 2 * self.matrixD[p, p]
            sumnotdiag += self.matrixS[p, i] * self.matrixD[p, p] * self.matrixS[p, j]

        return (sumdiag, sumnotdiag)

    def createSD(self):
        for i in range(3):
            for j in range(3):
                if i == j:
                    self.matrixD[i, j] = np.sign(self.matrixA[i, j] - self.sumforSD(i, j)[0])
                    self.matrixS[i, j] = np.sqrt(abs(self.matrixA[i, j] - self.sumforSD(i, j)[0]))  #math.sqrt no np.sqrt
                    print(f"d{i}{j} = {self.matrixD[i, j]} ")
                elif j > i:
                    self.matrixS[i, j] = (self.matrixA[i, j] - self.sumforSD(i, j)[1]) / (self.matrixD[i, i] * self.matrixS[i, i])
                print(f"s{i}{j} = {self.matrixS[ i , j ]} ")

    def _determinant_sqm(self)->float:
        det = 1
        for i in range(3):
            det *= self.matrixS[i, i] ** 2 * self.matrixD[i, i]
        return det

    def square_root_method(self):
        A = self.matrixA
        b = self.rightVector
        x = np.array([["x1"], ["x2"], ["x3"]])
        print(f"Starting system of equations A * x = b ")
        print("┌             ┐   ┌   ┐   ┌   ┐")
        print(f"│ {A[ 0 , 0 ]}   {A[ 0 , 1 ]}     {A[ 0 , 2 ]} │   │ {x[ 0 , 0 ]}│   │ {b[ 0 , 0 ]} │")
        print(f"│ {A[ 1 , 0 ]}   {A[ 1 , 1 ]}     {A[ 1 , 2 ]} │ * │ {x[ 1 , 0 ]}│ = │ {b[ 1 , 0 ]} │")
        print(f"│ {A[ 2 , 0 ]}    {A[ 2 , 1 ]}    {A[ 2 , 2 ]} │   │ {x[ 2 , 0 ]}│   │ {b[ 2 , 0 ]} │")
        print("└             ┘   └   ┘   └   ┘")

        print("\nCreating matrices S and D ")
        self.createSD()
        print(f"S = {self.matrixS} \n")
        print(f"D = {self.matrixD} \n")

        ST_mult_D = self.matrixS.T @ self.matrixD
        print(f"ST_mult_D = {ST_mult_D} \n")
        y = np.array([["y1"],["y2"],["y3"]])

        print(f"Finding y:  ST_D * y = b ")
        print(f"y1 = 3, y2 = 2, y3 = -2 => y_T = (3, 2, -2) ")
        print("┌                 ┐   ┌   ┐   ┌   ┐")
        print(f"│ {ST_mult_D[ 0 , 0 ]}   {ST_mult_D[ 0 , 1 ]}    {ST_mult_D[ 0 , 2 ]}│   │ {y[ 0 , 0 ]}│   │ {b[ 0 , 0 ]} │")
        print(f"│{ST_mult_D[ 1 , 0 ]}   {ST_mult_D[ 1 , 1 ]}    {ST_mult_D[ 1 , 2 ]}│ * │ {y[ 1 , 0 ]}│ = │ {b[ 1 , 0 ]} │")
        print(f"│{ST_mult_D[ 2 , 0 ]}    {ST_mult_D[ 2 , 1 ]}   {ST_mult_D[ 2 , 2 ]}│   │ {y[ 2 , 0 ]}│   │ {b[ 2 , 0 ]} │")
        print("└                 ┘   └   ┘   └   ┘\n")

        print("Final equation: S * x = y ")
        print("┌                   ┐   ┌   ┐   ┌   ┐")
        print(f"│ {self.matrixS[ 0 , 0 ]}   {self.matrixS[ 0 , 1 ]}     {self.matrixS[ 0 , 2 ]} │   │ {x[ 0 , 0 ]}│   │ {y[ 0 , 0 ]}│")
        print(f"│ {self.matrixS[ 1 , 0 ]}   {self.matrixS[ 1 , 1 ]}     {self.matrixS[ 1 , 2 ]} │ * │ {x[ 1 , 0 ]}│ = │ {y[ 1 , 0 ]}│")
        print(f"│ {self.matrixS[ 2 , 0 ]}    {self.matrixS[ 2 , 1 ]}    {self.matrixS[ 2 , 2 ]} │   │ {x[ 2 , 0 ]}│   │ {y[ 2 , 0 ]}│")
        print("└                   ┘   └   ┘   └   ┘")

        print(f"\n x1 = -1, x2= 4, x3 = -2 => x_T = (-1, 4, -2) \n")

        print(f"Determinant of matrix A = {self._determinant_sqm()} \n")
        print(f"Condition number = {self.condition_number()}")

    def check_minors(self, matrix)->bool:
        check = False
        det1 = matrix[0,0]
        print(f"Minor1 = {det1}")
        det2 = self.determinant(matrix[0 : 2 , 0 : 2])
        print(f"Minor2 = {det2}")
        det3 = self.determinant(matrix)
        print(f"Minor3 = {det3}")
        if det1 > 0 and det2 > 0 and det3 > 0:
            print("Second condition satisfied: A > 0")
            check = True
        return check



    def seidel_method(self):
        A = self.matrixA
        b = self.rightVector
        x = np.array([ [ "x1" ] , [ "x2" ] , [ "x3" ] ])
        x1 = x2 = x3 = 0
        eps = float(input("Enter eps (precision) : "))
        vecx_prev = np.array([[x1],[x2],[x3]]) #x0
        print(f"Starting approximation x0 = {vecx_prev} \n")
        vecx = np.array([[1],[1],[1]])
        val1 = float(abs(vecx[0,0] - vecx_prev[0,0]))
        val2 = float(abs(vecx[1,0] - vecx_prev[1,0]))
        val3 = float(abs(vecx[2,0] - vecx_prev[2,0]))
        inf_norm = max(val1, val2, val3)

        print(f"Starting system of equations A * x = b ")
        print("┌             ┐   ┌   ┐   ┌   ┐")
        print(f"│ {A[ 0 , 0 ]}   {A[ 0 , 1 ]}     {A[ 0 , 2 ]} │   │ {x[ 0 , 0 ]}│   │ {b[ 0 , 0 ]} │")
        print(f"│ {A[ 1 , 0 ]}   {A[ 1 , 1 ]}     {A[ 1 , 2 ]} │ * │ {x[ 1 , 0 ]}│ = │ {b[ 1 , 0 ]} │")
        print(f"│ {A[ 2 , 0 ]}    {A[ 2 , 1 ]}    {A[ 2 , 2 ]} │   │ {x[ 2 , 0 ]}│   │ {b[ 2 , 0 ]} │")
        print("└             ┘   └   ┘   └   ┘\n")

        print("Checking minors ")
        cond = self.check_minors(A)
        print(f"Second condition is {cond}, starting iterations \n")

        iterations = 0
        while inf_norm > eps:
            print(f"Iteration = { iterations }")

            vecx_prev = np.array([ [ x1 ] , [ x2 ] , [ x3 ] ])
            x1 = (-A[0,1] * x2 - A[0,2] * x3 + b[0,0]) / A[0,0]
            x2 = (-A[1,0] * x1 - A[1,2] * x3 + b[1,0]) / A[1,1]
            x3 = (-A[2,0] * x1 - A[2,1] * x2 + b[2,0]) / A[2,2]
            vecx = np.array([ [ x1 ] , [ x2 ] , [ x3 ] ])
            print(f"x1 = {x1}, x2 = {x2}, x3 = {x3}")

            print(f"checking if infinite norm of | x{iterations} - x{iterations - 1}| > eps ")
            val1 = float(abs(vecx[ 0 , 0 ] - vecx_prev[ 0 , 0 ]))
            val2 = float(abs(vecx[ 1 , 0 ] - vecx_prev[ 1 , 0 ]))
            val3 = float(abs(vecx[ 2 , 0 ] - vecx_prev[ 2 , 0 ]))
            inf_norm = max(val1 , val2 , val3)
            print(f"Infinite norm of | x{iterations} - x{iterations - 1}| = {inf_norm} \n")

            iterations += 1

        print(f"Solution = {vecx} \n with precision = {eps} \n ")
        print(f"Condition number = {self.condition_number()}")
        return vecx


    def determinant(self, matrix)->float:
        return np.linalg.det(matrix)


    def condition_number(self)->float:
        matrixA = self.matrixA
        A_inverse = np.linalg.inv(matrixA)
        print(f"inverse_matrix = {A_inverse} \n")

        A_abs = np.abs(matrixA)
        A_inverse_abs = np.abs(A_inverse)
        A_infnorm = max(float(np.sum(A_abs[0,0:])), float(np.sum(A_abs[1,0:])), float(np.sum(A_abs[2,0:])))
        A_inverse_infnorm = max(float(np.sum(A_inverse_abs[0,0:])), float(np.sum(A_inverse_abs[1,0:])), float(np.sum(A_inverse_abs[2,0:])))
        condition_num = A_infnorm * A_inverse_infnorm

        return condition_num

ex1 = Solution()
sol = ex1.seidel_method()
#sol = ex1.square_root_method()
#ex1.square_root_method()