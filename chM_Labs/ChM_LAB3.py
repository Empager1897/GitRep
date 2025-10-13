# Python 3.12
# Степеневий метод та обернений степеневий метод для симетричної додатньо визначеної матриці 3x3.
# Виводяться всі ітерації, проміжні обчислення, перевіряються умови застосовності та збіжності.

from math import sqrt
from typing import List, Tuple, Optional

Vector = List[float]
Matrix = List[List[float]]


def is_symmetric(A: Matrix, tol: float = 1e-12) -> bool:
    n = len(A)
    for i in range(n):
        if len(A[i]) != n:
            return False
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > tol:
                return False
    return True


def gershgorin_bounds(A: Matrix) -> Tuple[float, float]:
    # Дає грубі оцінки: усі власні значення у [min(c_i-R_i), max(c_i+R_i)]
    n = len(A)
    lowers = []
    uppers = []
    for i in range(n):
        c = A[i][i]
        R = sum(abs(A[i][j]) for j in range(n) if j != i)
        lowers.append(c - R)
        uppers.append(c + R)
    return min(lowers), max(uppers)


def cholesky_decompose(A: Matrix, tol: float = 1e-15) -> Matrix:
    # Розклад Холєцького A = L L^T; кидає виняток, якщо матриця не SPD
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                val = A[i][i] - s
                if val <= tol:
                    raise ValueError("Матриця не є додатньо визначеною (Cholesky провалився).")
                L[i][i] = sqrt(val)
            else:
                if abs(L[j][j]) <= tol:
                    raise ValueError("Ділення на нуль у Cholesky — матриця вироджена.")
                L[i][j] = (A[i][j] - s) / L[j][j]
    return L


def forward_substitution(L: Matrix, b: Vector) -> Vector:
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        s = sum(L[i][k] * y[k] for k in range(i))
        y[i] = (b[i] - s) / L[i][i]
    return y


def backward_substitution(U: Matrix, y: Vector) -> Vector:
    n = len(U)
    x = [0.0] * n
    for i in reversed(range(n)):
        s = sum(U[i][k] * x[k] for k in range(i + 1, n))
        x[i] = (y[i] - s) / U[i][i]
    return x


def solve_cholesky(L: Matrix, b: Vector) -> Vector:
    # Розв’язує A x = b, маючи L, де A = L L^T
    y = forward_substitution(L, b)
    x = backward_substitution(list(zip(*L)), y)  # L^T
    x = [float(v) for v in x]
    return x


def dot(u: Vector, v: Vector) -> float:
    return sum(ui * vi for ui, vi in zip(u, v))


def matvec(A: Matrix, x: Vector) -> Vector:
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def norm2(x: Vector) -> float:
    return sqrt(dot(x, x))


def normalize(x: Vector) -> Tuple[Vector, float]:
    nrm = norm2(x)
    if nrm == 0.0:
        raise ValueError("Норма вектора дорівнює нулю — нормування неможливе.")
    return [xi / nrm for xi in x], nrm


def rayleigh_quotient(A: Matrix, x: Vector) -> float:
    Ax = matvec(A, x)
    return dot(x, Ax) / dot(x, x)


def residual_norm(A: Matrix, x: Vector, lam: float) -> float:
    Ax = matvec(A, x)
    r = [Ax[i] - lam * x[i] for i in range(len(x))]
    return norm2(r)


def fmt_vec(x: Vector, prec: int = 10) -> str:
    return "[" + ", ".join(f"{xi:+.{prec}f}" for xi in x) + "]"


def power_method(
    A: Matrix,
    tol: float = 1e-12,
    max_iters: int = 100,
    x0: Optional[Vector] = None,
    verbose: bool = True,
) -> Tuple[float, Vector, int, bool]:
    n = len(A)
    if x0 is None:
        x = [1.0] * n
    else:
        if len(x0) != n:
            raise ValueError("Невірна розмірність початкового вектора.")
        x = x0[:]

    lam_prev: Optional[float] = None
    converged = False

    if verbose:
        print("\n=== Степеневий метод (максимальне власне значення) ===")
    for k in range(1, max_iters + 1):
        y = matvec(A, x)
        try:
            x_next, y_norm = normalize(y)
        except ValueError:
            raise ValueError("Отримано нульовий вектор на кроці нормування в степеневому методі.")

        lam = rayleigh_quotient(A, x_next)
        res = residual_norm(A, x_next, lam)
        delta = float("inf") if lam_prev is None else abs(lam - lam_prev)

        if verbose:
            print(f"Ітерація {k}:")
            print(f"  y = A x_k = {fmt_vec(y)}")
            print(f"  ||y||_2 = {y_norm:.12e}")
            print(f"  x_(k+1) = y / ||y|| = {fmt_vec(x_next)}")
            print(f"  λ_(k+1)  = {lam:.12e}")
            if lam_prev is not None:
                print(f"  |Δλ| = {delta:.12e}")
            print(f"  ||A x - λ x||_2 = {res:.12e}")

        # Умова зупинки: стабілізація λ і мала нев’язка
        if lam_prev is not None and delta < tol and res < tol:
            converged = True
            x = x_next
            lam_prev = lam
            if verbose:
                print("  -> Збіжність досягнута за заданим порогом.")
            break

        x = x_next
        lam_prev = lam

    if not converged and verbose:
        print("Попередження: досягнуто максимальну кількість ітерацій без досягнення порогу збіжності.")

    return lam_prev if lam_prev is not None else 0.0, x, k, converged


def inverse_power_method_min(
    A: Matrix,
    tol: float = 1e-12,
    max_iters: int = 100,
    x0: Optional[Vector] = None,
    verbose: bool = True,
) -> Tuple[float, Vector, int, bool]:
    # Для мінімального власного значення SPD-матриці використовуємо обернений степеневий метод через Cholesky.
    n = len(A)
    if x0 is None:
        x = [1.0] * n
    else:
        if len(x0) != n:
            raise ValueError("Невірна розмірність початкового вектора.")
        x = x0[:]

    # Cholesky (A = L L^T)
    L = cholesky_decompose(A)

    lam_prev: Optional[float] = None
    converged = False

    if verbose:
        print("\n=== Обернений степеневий метод (мінімальне власне значення) ===")
    for k in range(1, max_iters + 1):
        # Розв'язуємо A z = x (тобто z ≈ A^{-1} x)
        z = solve_cholesky(L, x)
        try:
            x_next, z_norm = normalize(z)
        except ValueError:
            raise ValueError("Отримано нульовий вектор на кроці нормування в оберненому степеневому методі.")

        lam = rayleigh_quotient(A, x_next)  # надійна оцінка λ_min для симетричної A
        res = residual_norm(A, x_next, lam)
        delta = float("inf") if lam_prev is None else abs(lam - lam_prev)

        if verbose:
            print(f"Ітерація {k}:")
            print(f"  Розв'язок A z = x_k -> z = {fmt_vec(z)}")
            print(f"  ||z||_2 = {z_norm:.12e}")
            print(f"  x_(k+1) = z / ||z|| = {fmt_vec(x_next)}")
            print(f"  λ_(k+1)  = {lam:.12e}")
            if lam_prev is not None:
                print(f"  |Δλ| = {delta:.12e}")
            print(f"  ||A x - λ x||_2 = {res:.12e}")

        if lam_prev is not None and delta < tol and res < tol:
            converged = True
            x = x_next
            lam_prev = lam
            """if verbose:
                print("  -> Збіжність досягнута за заданим порогом.")"""
            break

        x = x_next
        lam_prev = lam

    if not converged and verbose:
        print("Попередження: досягнуто максимальну кількість ітерацій без досягнення порогу збіжності.")

    return lam_prev if lam_prev is not None else 0.0, x, k, converged


def analyze_applicability(A: Matrix) -> None:
    print("=== Дослідження можливості застосування методів ===")
    print("Матриця A:")
    for row in A:
        print("  ", fmt_vec(row))

    sym = is_symmetric(A)
    print(f"\nПеревірка симетричності: {'ТАК' if sym else 'НІ'}")

    try:
        L = cholesky_decompose(A)
        spd = True
        print("Перевірка додатньої визначеності : ТАК")
    except Exception as e:
        spd = False
        print(f"Перевірка додатньої визначеності : НІ ({e})")

    low, up = gershgorin_bounds(A)
    print(f"Гершгорінові оцінки для спектра: всі λ(A) ∈ [{low:.6f}, {up:.6f}]")
    if low > 0:
        print("  Примітка: нижня гершгорінова межа > 0, тому всі власні значення гарантовано додатні (достатня умова).")
    else:
        print("  Примітка: гершгорінова нижня межа не гарантує додатність (умова лише достатня).")

    """print("\n Достатня умова збіжності степеневого методу:")
    print("- найбільше за модулем власне значення має бути простим;")
    print("- початковий вектор має ненульову проєкцію на власний вектор для найбільшого власного значення;")
    print("У програмі контролюється збіжність через |Δλ| і норму нев’язки ||A x - λ x||.")"""

    if not sym or not spd:
        print("\nУВАГА: методи надалі можуть бути незастосовними. Для пошуку λ_min обов’язково потрібна SPD.")
    print()


def read_matrix_from_stdin() -> Optional[Matrix]:
    try:
        print("Введіть матрицю 3×3 (три рядки по три числа, розділені пробілами).")
        rows: Matrix = []
        for i in range(3):
            line = input().strip()
            if not line:
                return None
            parts = line.replace(",", " ").split()
            if len(parts) != 3:
                return None
            rows.append([float(p) for p in parts])
        return rows
    except Exception:
        return None


def main():
    # Приклад SPD-матриці за замовчуванням
    A_default: Matrix = [
        [4.0, 1.0, 1.0],
        [1.0, 3.0, 0.0],
        [1.0, 0.0, 2.0],
    ]

    print("Обрати матрицю:")
    print("  1) Ввести вручну (3×3)")
    print("  2) Використати приклад за замовчуванням")
    choice = input("Введіть 1 або 2 і натисніть Enter: ").strip()

    if choice == "1":
        A = read_matrix_from_stdin()
        if A is None:
            print("Не вдалося зчитати матрицю. Використовується приклад за замовчуванням.")
            A = A_default
    else:
        A = A_default

    analyze_applicability(A)

    if not is_symmetric(A):
        print("Помилка: матриця не є симетричною. Степеневий метод для симетричної SPD матриці не застосовується.")
        return

    # Параметри методу
    tol = 1e-12
    max_iters = 100
    x0 = [1.0, 1.0, 1.0]  # початковий вектор

    # Максимальне власне значення (степеневий метод)
    lam_max, v_max, iters_max, conv_max = power_method(A, tol=tol, max_iters=max_iters, x0=x0, verbose=True)

    print("\n=== Результат для максимального власного значення ===")
    print(f"λ_max ≈ {lam_max:.12e}")
    print(f"Власний вектор (нормований): {fmt_vec(v_max)}")
    print(f"Ітерацій: {iters_max}, збіжність: {'ТАК' if conv_max else 'НІ'}")

    """# Мінімальне власне значення (обернений степеневий метод через Cholesky)
    try:
        lam_min, v_min, iters_min, conv_min = inverse_power_method_min(A, tol=tol, max_iters=max_iters, x0=x0, verbose=True)
        print("\n=== Результат для мінімального власного значення ===")
        print(f"λ_min ≈ {lam_min:.12e}")
        print(f"Власний вектор (нормований): {fmt_vec(v_min)}")
        print(f"Ітерацій: {iters_min}, збіжність: {'ТАК' if conv_min else 'НІ'}")
    except Exception as e:
        print("\nНе вдалося обчислити мінімальне власне значення оберненим степеневим методом:")
        print(f"  {e}")"""


if __name__ == "__main__":
    main()