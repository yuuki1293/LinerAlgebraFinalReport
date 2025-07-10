import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

# フォントエラーを回避するため、デフォルトフォントを使用
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

# データの読み込み
det_data = pd.read_csv('data/determinant_times.csv')
linear_data = pd.read_csv('data/linear_solver_times.csv')
eigen_data = pd.read_csv('data/eigenvalue_times.csv')

def power_law(x, a, b):
    """Power law: y = a * x^b"""
    return a * np.power(x, b)

def analyze_complexity(data, title, eng_title, eng_ylabel):
    """Analyze computational complexity and plot in English"""
    x = data.iloc[:, 0].values  # Size
    y = data.iloc[:, 1].values  # Time

    # Use only nonzero data
    mask = y > 0
    x_filtered = x[mask]
    y_filtered = y[mask]

    # Log transform
    log_x = np.log(x_filtered)
    log_y = np.log(y_filtered)

    # Linear regression for power law
    coeffs = np.polyfit(log_x, log_y, 1)
    power = coeffs[0]
    intercept = coeffs[1]
    a = np.exp(intercept)

    # Curve fitting
    popt, pcov = curve_fit(power_law, x_filtered, y_filtered, p0=[a, power])
    fitted_a, fitted_power = popt

    # R² value
    y_pred = power_law(x_filtered, fitted_a, fitted_power)
    ss_res = np.sum((y_filtered - y_pred) ** 2)
    ss_tot = np.sum((y_filtered - np.mean(y_filtered)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    print(f"\n{title} analysis result:")
    print(f"Estimated exponent: {fitted_power:.3f}")
    print(f"Coefficient a: {fitted_a:.6f}")
    print(f"R²: {r_squared:.6f}")

    # Theoretical values
    theoretical_powers = {
        "Determinant": 3.0,
        "Linear Solver": 3.0,
        "Eigenvalue": 3.0
    }
    theoretical = theoretical_powers.get(title, 3.0)
    print(f"Theoretical exponent: {theoretical}")
    print(f"Measured/Theoretical: {fitted_power/theoretical:.3f}")

    # Plot
    plt.figure(figsize=(12, 5))

    # Measured data
    plt.subplot(1, 2, 1)
    plt.scatter(x_filtered, y_filtered, alpha=0.6, label='Measured Data')
    x_fit = np.linspace(min(x_filtered), max(x_filtered), 100)
    plt.plot(x_fit, power_law(x_fit, fitted_a, fitted_power), 'r-',
             label=f'Fit: y = {fitted_a:.6f} * x^{fitted_power:.3f}')
    plt.xlabel('Matrix Size n')
    plt.ylabel(eng_ylabel)
    plt.title(f'{eng_title} - Computation Time')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Log-log plot
    plt.subplot(1, 2, 2)
    plt.scatter(log_x, log_y, alpha=0.6, label='Measured Data')
    log_x_fit = np.linspace(min(log_x), max(log_x), 100)
    plt.plot(log_x_fit, fitted_power * log_x_fit + np.log(fitted_a), 'r-',
             label=f'Fit: log(y) = {fitted_power:.3f} * log(x) + {np.log(fitted_a):.3f}')
    plt.xlabel('log(n)')
    plt.ylabel(f'log({eng_ylabel})')
    plt.title(f'{eng_title} - Log Plot')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'figures/{title}_complexity_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()  # Prevent memory leak

    return fitted_power, r_squared

# Analyze each computation
print("=== Complexity Analysis of Computation Time ===")

det_power, det_r2 = analyze_complexity(det_data, "Determinant", "Determinant", "Computation Time (ms)")
linear_power, linear_r2 = analyze_complexity(linear_data, "Linear Solver", "Linear Solver", "Computation Time (ms)")
eigen_power, eigen_r2 = analyze_complexity(eigen_data, "Eigenvalue", "Eigenvalue", "Computation Time (ms)")

# Summary
print("\n=== Summary of Analysis Results ===")
print(f"{'Computation':<15} {'Measured Exp':<12} {'Theory Exp':<12} {'Ratio':<8} {'R²':<8}")
print("-" * 60)
print(f"{'Determinant':<15} {det_power:<12.3f} {3.0:<12.1f} {det_power/3.0:<8.3f} {det_r2:<8.3f}")
print(f"{'Linear Solver':<15} {linear_power:<12.3f} {3.0:<12.1f} {linear_power/3.0:<8.3f} {linear_r2:<8.3f}")
print(f"{'Eigenvalue':<15} {eigen_power:<12.3f} {3.0:<12.1f} {eigen_power/3.0:<8.3f} {eigen_r2:<8.3f}")

# Detailed analysis
print("\n=== Detailed Analysis ===")
print("1. Determinant:")
print(f"   - Measured exponent: {det_power:.3f} (Theory: 3.0)")
print(f"   - Ratio: {det_power/3.0:.3f}")
print(f"   - Reason: Cache effect, compiler optimization")

print("\n2. Linear Solver:")
print(f"   - Measured exponent: {linear_power:.3f} (Theory: 3.0)")
print(f"   - Ratio: {linear_power/3.0:.3f}")
print(f"   - Reason: Same as determinant")

print("\n3. Eigenvalue:")
print(f"   - Measured exponent: {eigen_power:.3f} (Theory: 3.0)")
print(f"   - Ratio: {eigen_power/3.0:.3f}")
print(f"   - Reason: Iteration count of QR, overhead of convergence check")
