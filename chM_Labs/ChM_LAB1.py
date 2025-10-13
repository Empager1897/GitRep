"""import seaborn as sns
from sympy import symbols , Interval, diff, Derivative
from sympy.calculus.util import continuous_domain"""
import matplotlib.pyplot as plt
import numpy as np
import math

class Solution:
    def function(self, x):
        return x ** 2 + np.sin(x) - 12 * x - 0.25

    def phi_x1(self, x):
        return (x**2 + np.sin(x) - 0.25) / 12

    def phi_x1_der(self, x):
        return (2 * x + np.cos(x)) / 12

    def phi_x2(self , x) :
        argument = -x ** 2 + 12 * x + 0.25
        argument = np.clip(argument , -1 , 1)  # Keep arcsin argument in valid range
        return np.arcsin(argument)

    def phi_x2_der(self , x) :
        argument = -x ** 2 + 12 * x + 0.25

        # Handle the sqrt domain carefully
        sqrt_arg = 1 - argument ** 2
        # Use np.where to avoid invalid sqrt values
        denominator = np.where(sqrt_arg > 0 , np.sqrt(sqrt_arg) , np.nan)

        # Avoid division by zero or very small numbers
        result = (-2 * x + 12) / denominator
        # Replace extreme values with reasonable bounds
        result = np.clip(result , -1e10 , 1e10)  # Clip to avoid overflow
        return result

    def phi_x3(self , x) :
        argument = 12 * x - np.sin(x) + 0.25
        # Use np.where to handle negative arguments gracefully
        return np.where(argument >= 0 , np.sqrt(argument) , np.nan)

    def phi_x3_der(self , x) :
        argument = 12 * x - np.sin(x) + 0.25

        # Handle the case where argument is negative or zero
        sqrt_val = np.where(argument > 0 , np.sqrt(argument) , np.nan)

        # Avoid division by zero - use np.where for safe division
        denominator = np.where(np.abs(sqrt_val) > 1e-10 , 2 * sqrt_val , np.nan)

        result = np.where(denominator != 0 , (12 - np.cos(x)) / denominator , np.nan)
        return result

    def __init__(self , eps) :
        self._eps = eps
        self.functions = {
            "starting" : self.function ,
            "phi_x1" : self.phi_x1 ,
            "phi_x2" : self.phi_x2 ,
            "phi_x3" : self.phi_x3 ,
            "phi_x1_der" : self.phi_x1_der ,
            "phi_x2_der" : self.phi_x2_der ,
            "phi_x3_der" : self.phi_x3_der
        }

    def choose_interval_step(self)->tuple:
        a0 = float(input("Start: "))
        b0 = float(input("End: "))
        step = float(input("Step: "))
        return (a0, b0, step)

    def draw(self, function_name):
        print(f"Enter the parameters of the interval for x-axis for {function_name} function")
        sx , ex , stx = self.choose_interval_step()
        print(f"Enter the parameters of the interval for y-axis for {function_name} function")
        sy , ey , sty = self.choose_interval_step()
        x = np.arange(sx , ex , stx)
        y = self.functions[ function_name ](x)
        fig, axs = plt.subplots(1, 1, figsize=(12, 8), layout = 'constrained')
        axs.set_yticks(np.arange(sy, ey, sty))
        axs.set_xticks(np.arange(sx, ex, stx))
        axs.plot(x, y)
        axs.spines[ 'bottom' ].set_position('zero')
        axs.spines[ 'left' ].set_position('zero')
        axs.spines[ 'right' ].set_visible(False)
        axs.spines[ 'top' ].set_visible(False)
        axs.xaxis.set_ticks_position('bottom')
        axs.yaxis.set_ticks_position('left')
        plt.show()

    def dichotomy(self):
        self.draw("starting")
        print("Choose the interval for the root of the equation ")
        a, b, step = self.choose_interval_step()  # -1, 0.5, 1 //  12, 12.1, 1
        if self.function(a) * self.function(b) >= 0:
            return "Wrong interval"
        step_counter = 0
        x = (a + b) / 2.0 #x0

        if input("Is function continuos?(t/f) ") == "t":
            print(f"\nApriori assessment of number of iterations >= {int(math.log2((b - a) / self._eps))} \n")
            while True:
                x_prev = x
                f_x_prev = self.function(x_prev)
                f_a = self.function(a)
                f_b = self.function(b)

                if self.function(a) * self.function(x_prev) < 0:
                    b = x # f(a)
                elif self.function(b) * self.function(x_prev) < 0:
                    a = x # f(b)
                if self.function(x) == 0:
                    print("x - solution")
                    break

                x = (a + b) / 2.0 # f(x), xn
                f_x = self.function(x)
                step_counter += 1

                print(f"Current iteration: {step_counter}")
                print(f"Current interval (a,b) = {(a , b)}")
                print(f"x{step_counter} = {x}")
                print(f"x{step_counter - 1} = {x_prev}")
                print(f"f(a) = {f_a}")
                print(f"f(b) = {f_b}")
                print(f"f(x{step_counter - 1}) = {f_x_prev}")
                print(f"f(x{step_counter}) = {f_x}")

                if abs(x - x_prev) <= self._eps:
                    print("Stop condition reached")
                    break

                print("\n")

            print(f"Aposteriori assessment of number of iterations = {step_counter} ")
        else:
            print("Function is not continuos or wrong interval (f(a) * f(b) >= 0)")
        return x

    def simple_iteration(self):
        self.draw("starting")
        print("Choose the interval for the root of the equation ")
        a, b, step = self.choose_interval_step() # -0.5 0.5 0.1
        step_counter = 0
        x = (a + b) / 2.0
        x_prev = x # x0
        phi = None
        beta = max(abs(a - x_prev), abs(b - x_prev))
        q = None #max(abs(phi_x_der)) < 1
        funcs = ["phi_x1", "phi_x2", "phi_x3"]
        derivs = ["phi_x1_der", "phi_x2_der", "phi_x3_der"]
                                                # 0 13 0.5 / 0 10 1
        for func, deriv in zip(funcs, derivs):
            self.draw(func)
            self.draw(deriv)
            check = input("Are the function and it`s derivative suitable? ")
            if check == "ok":
                q = abs(float(input(f"Enter the max absolute value (q) of the derivative for {deriv} "))) #0.1
                if q < 1 and abs(self.functions[func](x_prev) - x_prev) < (1 - q) * beta:
                    phi = func
                    break

        print(f"\n Apriori assessment of number of iterations >= {int(math.log(abs(self.functions[phi](x_prev) - x_prev)/
                                                (self._eps - self._eps * q)) / math.log(1 / q)) + 1}\n ")
        while q < 1 :
            step_counter += 1
            x_prev = x
            x = self.functions[phi](x_prev)

            print(f"Current iteration: {step_counter}")
            print(f"x{step_counter - 1} = {x_prev}")
            print(f"x{step_counter} = {x}")
            print("\n")

            if q < 0.5 and abs(x - x_prev) <= (1-q) * self._eps / q:
                break
            if q >= 0.5 and abs(x - x_prev) <= self._eps:
                 break

        print(f"Aposteriori assessment of number of iterations = {step_counter}\n")
        return x

# Solutions 12.06085205078125 ; -0.02264404296875
ex1 = Solution(0.0001)
sol = ex1.simple_iteration()
#sol = ex1.dichotomy()
print(f"Solution: {sol}")
print(f"f(sol) = {ex1.function(sol)}")
#ex1.draw("starting")

"""Start: -0.5
End: 0.5
Step: 0.1
Enter the parameters of the interval for x-axis for phi_x1 function
Start: -1
End: 1
Step: 0.05
Enter the parameters of the interval for y-axis for phi_x1 function
Start: -0.5
End: 0.5
Step: 0.05
Enter the parameters of the interval for x-axis for phi_x1_der function
Start: -1
End: 1
Step: 0.05
Enter the parameters of the interval for y-axis for phi_x1_der function
Start: -0.5
End: 0.5
Step: 0.05"""